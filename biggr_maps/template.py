import json
import math
from typing import TextIO
from biggr_maps import map


def load_as_template(fp: TextIO) -> map.Map:
    data = json.load(fp)
    m = map.Map(
        name=data[0]["map_name"],
        description=data[0]["map_description"],
        homepage=data[0].get("homepage"),
    )

    nodes = {k: map.node_from_dict(v) for k, v in data[1]["nodes"].items()}

    reactions = []
    for reaction_data in data[1]["reactions"].values():
        name = reaction_data["name"]
        bigg_id = reaction_data["bigg_id"]
        reversibility = reaction_data.get("reversibility", True)
        label_x = reaction_data["label_x"]
        label_y = reaction_data["label_y"]

        angle = 0

        associated_metabolites = {}
        mid_marker = None
        plus_multi_marker = None
        minus_multi_marker = None
        for segment_id, segment in reaction_data.get("segments", {}).items():
            for k in ["from_node_id", "to_node_id"]:
                node = nodes[segment[k]]
                if node.node_type == "metabolite":
                    if k == "from_node_id":
                        b1 = segment.get("b2")
                        b2 = segment.get("b1")
                    else:
                        b1 = segment.get("b1")
                        b2 = segment.get("b2")
                    if b1 is not None:
                        b1 = (b1["x"], b1["y"])
                    if b2 is not None:
                        b2 = (b2["x"], b2["y"])
                    associated_metabolites[node.bigg_id] = (node, (b1, b2))
                    if node.node_is_primary:
                        if k == "from_node_id":
                            other_node = nodes[segment["to_node_id"]]
                        else:
                            other_node = nodes[segment["from_node_id"]]
                        coefficient = next(
                            x["coefficient"]
                            for x in reaction_data["metabolites"]
                            if x["bigg_id"] == node.bigg_id
                        )
                        angle = math.atan2(node.y - other_node.y, node.x - other_node.x)
                        if coefficient < 0:
                            angle = math.remainder(angle + math.pi, 2 * math.pi)
                        
                        if other_node.node_type == "multimarker" or other_node.node_type == "midmarker":
                            if coefficient < 0:
                                minus_multi_marker = other_node
                            else:
                                plus_multi_marker = other_node
                elif node.node_type == "midmarker":
                    mid_marker = node
        if mid_marker is None:
            print("Requires one mid marker.")
            continue

        reaction = map.AutoReactionWithOptionalMetabolites(
            name=name,
            bigg_id=bigg_id,
            mid_marker=mid_marker,
            angle=angle,
            reversibility=reversibility,
            label_x=label_x,
            label_y=label_y,
            minus_multi_marker=minus_multi_marker,
            plus_multi_marker=plus_multi_marker,
        )
        for metabolite_data in reaction_data["metabolites"]:
            node, b1_b2 = associated_metabolites[metabolite_data["bigg_id"]]
            if node.node_is_primary:
                reaction.add_metabolite(
                    node=node,
                    coefficient=metabolite_data["coefficient"],
                    b1_b2=b1_b2
                )
            else:
                reaction.add_optional_metabolite(node=node, coefficient=coefficient, b1_b2=b1_b2)
        reactions.append(reaction)
    
    for label_data in data[1]["text_labels"].values():
        label = map.TextLabel(**label_data)
        m.add_label(label)
    return m, reactions, nodes
