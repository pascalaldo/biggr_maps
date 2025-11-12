import math
from biggr_maps import map, pathway
from itertools import islice, pairwise
import xml.etree.ElementTree as ET
from functools import partial


# Not available in python 3.10 yet, should be removed and
# itertools.batched should be used.
def batched(iterable, n, *, strict=False):
    # batched('ABCDEFG', 2) → AB CD EF G
    if n < 1:
        raise ValueError("n must be at least one")
    iterator = iter(iterable)
    while batch := tuple(islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError("batched(): incomplete batch")
        yield batch

def is_within(coord, refcoord1, refcoord2):
    min_x = refcoord1[0] if refcoord1[0] < refcoord2[0] else refcoord2[0]
    min_y = refcoord1[1] if refcoord1[1] < refcoord2[1] else refcoord2[1]
    max_x = refcoord1[0] if refcoord1[0] >= refcoord2[0] else refcoord2[0]
    max_y = refcoord1[1] if refcoord1[1] >= refcoord2[1] else refcoord2[1]

    if min_x == max_x:
        return (coord[0] == min_x) and (min_y <= coord[1] <= max_y)
    elif min_y == max_y:
        return (coord[1] == min_y) and (min_x <= coord[1] <= max_x)
    return False


def kgml_to_escher_map(filename: str, scale_factor: float=3.0) -> map.Map:
    tree = ET.parse(filename)
    root = tree.getroot()

    name = root.attrib["title"]
    m = map.Map(name, description=name)

    nodes = {}
    reactions = {}
    reaction_synonyms = {}

    for child in root:
        if child.tag != "entry":
            continue
        if child.attrib.get("type") != "compound":
            continue
        child_id = int(child.attrib["id"])
        node_name = child.attrib["name"].removeprefix("cpd:")
        graphic = None
        for entry_child in child:
            if entry_child.tag != "graphics":
                continue
            if entry_child.attrib.get("type") != "circle":
                continue
            graphic = entry_child
            break
        if graphic is None:
            continue
        x = float(graphic.attrib["x"]) * scale_factor
        y = float(graphic.attrib["y"]) * scale_factor
        metabolite_node = map.MetaboliteNode(
            bigg_id=node_name, name=node_name, x=x, y=y, node_is_primary=True
        )
        nodes[child_id] = metabolite_node

    reaction_lines = {}
    for child in root:
        if child.tag != "entry":
            continue
        if child.attrib.get("type") not in ["gene", "ortholog", "reaction"]:
            continue
        reaction_names = [
            x.removeprefix("rn:") for x in child.attrib["reaction"].split(" ")
        ]
        reaction_name = reaction_names[0]
        for rn in reaction_names:
            reaction_synonyms[rn] = reaction_name
        lines = []
        for entry_child in child:
            if entry_child.tag != "graphics":
                continue
            if entry_child.attrib.get("type") != "line":
                continue
            coords = entry_child.attrib["coords"]
            coords = list(batched((float(x) * scale_factor for x in coords.split(",")), 2))
            lines.append(coords)
        reaction_lines[reaction_name] = lines

    for child in root:
        if child.tag != "reaction":
            continue
        reaction_id = int(child.attrib["id"])
        reaction_name = (
            child.attrib["name"].split(" ", maxsplit=1)[0].removeprefix("rn:")
        )
        reaction_name = reaction_synonyms[reaction_name]
        if reactions.get(reaction_name) is not None:
            continue
        metabolites = []
        for entry_child in child:
            if entry_child.tag == "substrate":
                coefficient = -1
            elif entry_child.tag == "product":
                coefficient = 1
            else:
                continue
            node = nodes[int(entry_child.attrib["id"])]
            metabolites.append({"coefficient": coefficient, "node": node})

        lines = reaction_lines[reaction_name]
        if lines:
            overlapping_line = lines[0]
            flip = False
            for metabolite_info in metabolites:
                node = metabolite_info["node"]
                if ((node.x - overlapping_line[0][0])**2 + (node.y - overlapping_line[0][1])**2) <= (25 * scale_factor)**2:
                    flip = metabolite_info["coefficient"] > 0
                    print("Found flip")
                    break
            for line in lines[1:]:
                overlapping_segments = []
                for coord0_1, coord0_2 in pairwise(line):
                    for coord1_1, coord1_2 in pairwise(overlapping_line):
                        if is_within(coord0_1, coord1_1, coord1_2):
                            if is_within(coord0_2, coord1_1, coord1_2):
                                overlapping_segments.append((coord0_1, coord0_2))
                            elif is_within(coord1_1, coord0_1, coord0_2):
                                overlapping_segments.append((coord0_1, coord1_1))
                            elif is_within(coord1_2, coord0_1, coord0_2):
                                overlapping_segments.append((coord0_1, coord1_2))
                        elif is_within(coord0_2, coord1_1, coord1_2):
                            if is_within(coord1_1, coord0_1, coord0_2):
                                overlapping_segments.apend((coord0_2, coord1_1))
                            elif is_within(coord1_2, coord0_1, coord0_2):
                                overlapping_segments.apend((coord0_2, coord1_2))
                ol = []
                for segment in overlapping_segments:
                    for o in ol:
                        if o[0][0] == segment[0][0] and o[0][1] == segment[0][1]:
                            if (o[0][0] == o[1][0] == segment[1][0]) or (o[0][1] == o[1][1] == segment[1][1]):
                                o.pop(0)
                            o.insert(0, segment[1])
                            break
                        elif o[0][0] == segment[1][0] and o[0][1] == segment[1][1]:
                            if (o[0][0] == o[1][0] == segment[0][0]) or (o[0][1] == o[1][1] == segment[0][1]):
                                o.pop(0)
                            o.insert(0, segment[0])
                            break
                        elif o[1][0] == segment[0][0] and o[1][1] == segment[0][1]:
                            if (o[0][0] == o[1][0] == segment[1][0]) or (o[0][1] == o[1][1] == segment[1][1]):
                                o.pop(-1)
                            o.append(segment[1])
                            break
                        elif o[1][0] == segment[1][0] and o[1][1] == segment[1][1]:
                            if (o[0][0] == o[1][0] == segment[0][0]) or (o[0][1] == o[1][1] == segment[0][1]):
                                o.pop(-1)
                            o.append(segment[0])
                            break
                    else:
                        ol.append([segment[0], segment[1]])
                if not ol:
                    print(f"No overlapping line found for reaction {reaction_name}")
                    overlapping_line = lines[0]
                    break
                else:
                    overlapping_line = ol[0]
            lengths = [
                ((coord2[0] - coord1[0]) ** 2 + (coord2[1] - coord1[1]))
                for coord1, coord2 in pairwise(overlapping_line)
            ]
            index_max = max(range(len(lengths)), key=lengths.__getitem__)
            dx = (overlapping_line[index_max + 1][0] - overlapping_line[index_max][0])
            dy = (overlapping_line[index_max + 1][1] - overlapping_line[index_max][1])
            mid_x = overlapping_line[index_max][0] + dx/2
            mid_y = overlapping_line[index_max][1] + dy/2     
            angle = math.atan2(dy, dx)
            if flip:
                angle = math.remainder(angle + math.pi, 2 * math.pi)
            mid_marker = map.MidMarkerNode(mid_x, mid_y)
            
            reaction = map.AutoReactionWithOptionalMetabolites(
                bigg_id=reaction_name,
                name=reaction_name,
                mid_marker=mid_marker,
                angle=angle,
                unit=50,
            )
            for metabolite_info in metabolites:
                reaction.add_metabolite(**metabolite_info)
            reactions[reaction_name] = reaction
        else:
            for metabolite_info in metabolites:
                m.add_node(metabolite_info["node"])
            mid_markers = [r.mid_marker for r in reactions.values()]
            reaction = pathway.place_reaction_on_backbone(
                map=m,
                name=reaction_name,
                bigg_id=reaction_name,
                reaction_info=[(x["coefficient"], x["node"]) for x in metabolites],
                additional_mid_markers=mid_markers,
                placement_f=partial(pathway.alternating_pathways_sides, spacing=100)
            )
            reactions[reaction_name] = reaction
    return m, reactions, nodes