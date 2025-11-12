from typing import Any, Dict, List, Optional, Tuple
from biggr_maps.map import Map, AutoReaction, MetaboliteNode, MidMarkerNode, Node
import math


def alternating_pathways_sides(
    map: Map,
    node1: Node,
    node2: Node,
    mid_markers: List[MidMarkerNode],
    spacing: float = 200,
    centered: bool = True,
):
    dx = node2.x - node1.x
    dy = node2.y - node1.y
    distance = math.sqrt(dx**2 + dy**2)
    center_point = (node1.x + dx/2, node1.y + dy/2)
    normal = (-dy/distance, dx/distance)

    i = 1 if centered else 0
    while True:
        side = i%2
        f = (i//2)
        if not centered:
            f += 0.5
        mid_marker_xy = (
            center_point[0] + (1 if side else -1) * f * spacing * normal[0],
            center_point[1] + (1 if side else -1) * f * spacing * normal[1],
        )
        if not any(math.sqrt((mid_marker_xy[0] - node.x)**2 + (mid_marker_xy[1] - node.y)**2) < spacing*0.99 for node in mid_markers):
            return mid_marker_xy
        i += 1

def place_reaction_on_backbone(
    map: Map,
    name: str,
    bigg_id: str,
    reaction_info: List[Tuple[float, MetaboliteNode]],
    placement_f = alternating_pathways_sides,
    add_metabolite_opts = None,
    additional_mid_markers: Optional[List[MidMarkerNode]] = None,
    **kwargs
):
    if add_metabolite_opts is None:
        add_metabolite_opts = {}
    backbone_nodes = [
        (coeff, n) for coeff, n in reaction_info if n.identifier in map.nodes
    ]
    while len(backbone_nodes) > 1:
        if (backbone_nodes[0][0] > 0) == (backbone_nodes[1][0] > 0):
            backbone_nodes.pop(1)
        else:
            break
    if len(backbone_nodes) < 2:
        raise ValueError(
            f"Two metabolite nodes should be present in the map already, found {len(backbone_nodes)}."
        )

    if (backbone_nodes[0][0] > 0) == (backbone_nodes[1][0] > 0):
        raise ValueError(f"The two nodes appear at the same side of the reaction.")

    if backbone_nodes[0][0] > 0:
        plus_node = backbone_nodes[0][1]
        minus_node = backbone_nodes[1][1]
    else:
        plus_node = backbone_nodes[1][1]
        minus_node = backbone_nodes[0][1]

    angle = math.atan2(plus_node.y - minus_node.y, plus_node.x - minus_node.x)
    # mid_marker = MidMarkerNode(
    #     plus_node.x + (minus_node.x - plus_node.x) / 2,
    #     plus_node.y + (minus_node.y - plus_node.y) / 2,
    # )
    mid_markers = [node for node in map.nodes.values() if node.node_type == "midmarker"]
    if additional_mid_markers is not None:
        mid_markers.extend(additional_mid_markers)
    mid_marker = MidMarkerNode(*placement_f(map, plus_node, minus_node, mid_markers))

    reaction = AutoReaction(
        name=name, bigg_id=bigg_id, mid_marker=mid_marker, angle=angle, **kwargs
    )
    for coeff, met in reaction_info:
        reaction.add_metabolite(met, coeff, **add_metabolite_opts)

    return reaction
