import networkit as nk
import momepy
import osmnx as ox
import networkx as nx

from osmnx import settings


def load_and_save_graph(city: str, country: str, instance_dir: str):
    location = f"{city}, {country}"
    # keep smoothness property of the road (osm tag)
    ox.settings.useful_tags_way = settings.useful_tags_way + ["smoothness"]
    # download location with OSMnx
    G = ox.graph_from_place(
        location,
        network_type="drive",
        custom_filter='["highway"]["highway"!~"footway|bridleway|path|service|track|steps"]',
    )

    G.remove_edges_from(nx.selfloop_edges(G))

    # graph to gdf
    gdf = ox.graph_to_gdfs(G, nodes=False, edges=True)

    # map nan to empty best value
    gdf["smoothness"] = gdf["smoothness"].fillna("excellent")

    smoothness = [
        "excellent",
        "good",
        "intermediate",
        "bad",
        "very_bad",
        "horrible",
        "very_horrible",
        "impassable",
    ]
    # some edges have an array of smoothness values "['bad', 'good']",
    # map them to the average of the values
    gdf["smoothness"] = gdf["smoothness"].apply(
        lambda x: sum([0 if i not in smoothness else smoothness.index(i) for i in x])
        / len(x)
        if isinstance(x, list)
        else 0
        if x not in smoothness
        else smoothness.index(x)
    )

    # clean up dataframe
    gdf = gdf.rename(
        columns={"osmid": "wayid", "geometry": "geometry", "smoothness": "health"}
    )
    gdf = gdf.reset_index()
    gdf = gdf[["wayid", "geometry", "health", "u", "v"]]
    gdf.wayid = gdf.wayid.apply(lambda x: x if type(x) == list else [x])
    gdf = gdf.to_crs(3857)

    gdf.to_csv("tkp.csv")

    # convert to networkx graph using momepy and set name of length
    nx_graph = momepy.gdf_to_nx(gdf, length="length", multigraph=False)

    nk_graph = nk.nxadapter.nx2nk(
        nx_graph, data=True, typeMap={"wayid": str, "u": str, "v": str}
    )

    # edges = momepy.nx_to_gdf(nx_graph, points=False, lines=True)

    # nk_graph.indexEdges()
    # state = nk_graph.attachEdgeAttribute("state", float)
    # cost = nk_graph.attachEdgeAttribute("cost", float)
    # wayids = nk_graph.attachEdgeAttribute("wayids", str)
    # u = nk_graph.attachEdgeAttribute("sourceid", str)
    # v = nk_graph.attachEdgeAttribute("destid", str)

    # def setEdgeAttributes(row):
    #     # print(row.name)
    #     state[row.name] = row.health
    #     cost[row.name] = row.length
    #     wayids[row.name] = str(row.wayid)
    #     u[row.name] = str(row.u)
    #     v[row.name] = str(row.v)

    # edges.apply(setEdgeAttributes, axis=1)

    state = nk_graph.getEdgeAttribute("health", float)
    cost = nk_graph.getEdgeAttribute("length", float)
    wayids = nk_graph.getEdgeAttribute("wayid", str)
    u = nk_graph.getEdgeAttribute("u", str)
    v = nk_graph.getEdgeAttribute("v", str)

    instance_name = location.replace(",", "-").replace(" ", "")
    state.write(f"{instance_dir}/{instance_name}.state")
    cost.write(f"{instance_dir}/{instance_name}.cost")
    wayids.write(f"{instance_dir}/{instance_name}.wayids")
    u.write(f"{instance_dir}/{instance_name}.sourceid")
    v.write(f"{instance_dir}/{instance_name}.destid")
    nk.writeGraph(
        nk_graph, f"{instance_dir}/{instance_name}.nkb", nk.Format.NetworkitBinary
    )


if __name__ == "__main__":
    instance_dir = "instances"

    load_and_save_graph("Berlin", "Germany", instance_dir)
    load_and_save_graph("Mitte", "Berlin, Germany", instance_dir)
    load_and_save_graph("Treptow-Köpenick", "Berlin, Germany", instance_dir)
