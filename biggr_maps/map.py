from dataclasses import dataclass
import math
from typing import Any, Dict, List, Optional, Tuple, Union


class Map:
    def __init__(
        self,
        name: str,
        description: str,
        homepage: Optional[str] = None,
        schema: Optional[str] = None,
        canvas: Optional[Tuple[float, float, float, float]] = None,
    ):
        if homepage is None:
            homepage = "https://escher.github.io"
        if schema is None:
            schema = "https://escher.github.io/escher/jsonschema/1-0-0#"
        if canvas is None:
            canvas = (0, 0, 1000, 1000)
        self.name = name
        self.description = description
        self.homepage = homepage
        self.schema = schema
        self.reactions = {}
        self.nodes = {}
        self.segments = {}
        self.text_labels = {}
        self.canvas = canvas
        self._node_counter = 0
        self._reaction_counter = 0
        self._segment_counter = 0

    def add_node(self, node: Optional["Node"]):
        if node is None or node.identifier is not None:
            return
        while self._node_counter in self.nodes:
            self._node_counter += 1
        node.identifier = self._node_counter
        self.nodes[node.identifier] = node
        self._node_counter += 1
    
    def add_segment(self, segment: "Segment"):
        while self._segment_counter in self.segments:
            self._segment_counter += 1
        segment.identifier = self._segment_counter
        self.segments[segment.identifier] = segment
        self._segment_counter += 1

    def add_reaction(self, reaction: "Reaction"):
        while self._reaction_counter in self.reactions:
            self._reaction_counter += 1
        reaction.identifier = self._reaction_counter
        reaction._map = self
        self.reactions[reaction.identifier] = reaction
        self._reaction_counter += 1

        for _, node in reaction.metabolites:
            self.add_node(node)
        self.add_node(reaction.mid_marker)
        self.add_node(reaction.multi_markers[0])
        self.add_node(reaction.multi_markers[1])

        for segment in reaction.segments:
            self.add_segment(segment)

    def fit_canvas(self, spacing: float = 100, expand_only=False):
        if expand_only:
            min_x = self.canvas[0]
            min_y = self.canvas[1]
            max_x = self.canvas[0] + self.canvas[2]
            max_y = self.canvas[1] + self.canvas[3]
        else:
            min_x, max_x, min_y, max_y = 0, 0, 0, 0
        for node in self.nodes.values():
            min_x = min(min_x, node.x)
            max_x = max(max_x, node.x)
            min_y = min(min_y, node.y)
            max_y = max(max_y, node.y)
        min_x -= spacing
        min_y -= spacing
        max_x += spacing
        max_y += spacing
        self.canvas = (min_x, min_y, max_x - min_x, max_y - min_y)

    def to_escher(self):
        d_header = {
            "map_name": self.name,
            "map_description": self.description,
            "homepage": self.homepage,
            "schema": self.schema,
        }
        d_body = {
            "reactions": {k: v.to_escher() for k, v in self.reactions.items()},
            "nodes": {k: v.to_escher() for k, v in self.nodes.items()},
            "text_labels": {k: v.to_escher() for k, v in self.text_labels.items()},
            "canvas": {
                "x": self.canvas[0],
                "y": self.canvas[1],
                "width": self.canvas[2],
                "height": self.canvas[3],
            },
        }
        return [d_header, d_body]


class Node:
    def __init__(self, x: Optional[float] = None, y: Optional[float] = None):
        self.identifier = None
        self.node_type = None
        self.x = x
        self.y = y

    def to_escher(self):
        d = self.__dict__.copy()
        del d["identifier"]
        return d

    def copy(self):
        o = self.__class__(**{k: v for k, v in self.__dict__.items() if k not in ["node_type", "identifier"]})
        o.node_type = self.node_type
        return o
    
class MetaboliteNode(Node):
    def __init__(
        self,
        bigg_id: str,
        name: str,
        x: Optional[float] = None,
        y: Optional[float] = None,
        label_x: Optional[float] = None,
        label_y: Optional[float] = None,
        node_is_primary: bool = False,
    ):
        super().__init__(x, y)
        self.node_type = "metabolite"
        self.bigg_id = bigg_id
        self.name = name
        self.label_x = label_x
        self.label_y = label_y
        self.node_is_primary = node_is_primary


class MultiMarkerNode(Node):
    def __init__(self, x: float, y: float):
        super().__init__(x, y)
        self.node_type = "multimarker"


class MidMarkerNode(Node):
    def __init__(self, x: float, y: float):
        super().__init__(x, y)
        self.node_type = "midmarker"


class Segment:
    def __init__(
        self,
        from_node: Node,
        to_node: Node,
        b1: Optional[Tuple[float, float]] = None,
        b2: Optional[Tuple[float, float]] = None,
    ):
        self.identifier = None
        self.from_node = from_node
        self.to_node = to_node
        self.b1 = b1
        self.b2 = b2

    def to_escher(self):
        d = {
            "from_node_id": self.from_node.identifier,
            "to_node_id": self.to_node.identifier,
            "b1": None if self.b1 is None else {"x": self.b1[0], "y": self.b1[1]},
            "b2": None if self.b2 is None else {"x": self.b2[0], "y": self.b2[1]},
        }
        return d

def node_from_dict(d: Dict[str, Any]) -> Node:
    node_class = Node
    if d["node_type"] == "metabolite":
        node_class = MetaboliteNode
    elif d["node_type"] == "midmarker":
        node_class = MidMarkerNode
    elif d["node_type"] == "multimarker":
        node_class = MultiMarkerNode
    return node_class(
        **{k: v for k, v in d.items() if k != "node_type"}
    )

class Reaction:
    def __init__(
        self,
        name: str,
        bigg_id: str,
        label_x: float,
        label_y: float,
        mid_marker: MidMarkerNode,
        plus_multi_marker: Optional[MultiMarkerNode],
        minus_multi_marker: Optional[MultiMarkerNode],
        reversibility: bool = True,
        gene_reaction_rule: Optional[str] = None,
        genes: Optional[List[Dict[str, str]]] = None,
    ):
        self.identifier = None
        self._map: Optional[Map] = None
        self.name = name
        self.bigg_id = bigg_id
        self.label_x = label_x
        self.label_y = label_y
        self.reversibility = reversibility
        self.mid_marker = mid_marker
        self.multi_markers = (minus_multi_marker, plus_multi_marker)
        self.metabolites = []
        self.segments = []
        self.gene_reaction_rule = gene_reaction_rule
        self.genes = genes
        self._add_multi_marker_segments()

    def add_segment(self, segment: "Segment"):
        if self._map is not None:
            self._map.add_segment(segment)
        self.segments.append(segment)
    
    def _add_multi_marker_segments(self):
        if self.multi_markers[0] is not None:
            self.add_segment(Segment(self.multi_markers[0], self.mid_marker))
        if self.multi_markers[1] is not None:
            self.add_segment(Segment(self.mid_marker, self.multi_markers[1]))

    def add_metabolite(self, node: MetaboliteNode, coefficient: Union[float, int]):
        self.metabolites.append((coefficient, node))
        if self._map is not None:
            self._map.add_node(node)

        ref_node = self.mid_marker
        plus_minus = coefficient > 0
        if self.multi_markers[plus_minus] is not None:
            ref_node = self.multi_markers[plus_minus]

        if plus_minus:
            segment = Segment(ref_node, node)
        else:
            segment = Segment(node, ref_node)
        self.add_segment(segment)

    def to_escher(self):
        d = {
            "name": self.name,
            "bigg_id": self.bigg_id,
            "reversibility": self.reversibility,
            "label_x": self.label_x,
            "label_y": self.label_y,
            "metabolites": [
                {"coeficient": coeff, "bigg_id": node.bigg_id}
                for coeff, node in self.metabolites
            ],
            "segments": {
                segment.identifier: segment.to_escher() for segment in self.segments
            },
        }
        if self.gene_reaction_rule is not None:
            d["gene_reaction_rule"] = self.gene_reaction_rule
        if self.genes is not None:
            d["genes"] = self.genes
        return d


def cubic_bezier_bt(t, b0, b1, b2, b3):
    bt = (
        ((1 - t) ** 3) * b0
        + 3 * t * ((1 - t) ** 2) * b1
        + 3 * (t**2) * (1 - t) * b2
        + (t**3) * b3
    )
    return bt


def non_primary_scaling(x):
    return (1 - (min(x - 1, 5) / 5)) * 0.3 + 0.5


def default_text_offset(x):
    return 20 + x * 12


class PlacementOptions:
    def __init__(
        self,
        delta=math.pi * 0.15,
        delta_tolerance=0.5,
        no_primary_length_f=None,
        scale=3.0,
        b1_scale=0.3,
        b2_scale=0.8,
        text_y_correction=6,
        text_offset_f=None,
        placement_f=None,
    ):
        self.delta = delta
        self.delta_tolerance = delta_tolerance
        self.scale = scale
        self.b1_scale = b1_scale
        self.b2_scale = b2_scale
        self.text_y_correction = text_y_correction

        self.placement_f = placement_f
        if no_primary_length_f is None:
            self.no_primary_length_f = non_primary_scaling
        else:
            self.no_primary_length_f = no_primary_length_f
        if text_offset_f is None:
            self.text_offset_f = default_text_offset
        else:
            self.text_offset_f = text_offset_f


class AutoReaction(Reaction):
    def __init__(
        self,
        bigg_id: str,
        mid_marker: MidMarkerNode,
        angle: float,
        unit: float = 50,
        text_y_correction: float = 8,
        label_x: Optional[float] = None,
        label_y: Optional[float] = None,
        minus_multi_marker: Optional[Union[MultiMarkerNode, MidMarkerNode]] = None,
        plus_multi_marker: Optional[Union[MultiMarkerNode, MidMarkerNode]] = None,
        **kwargs
    ):
        self.angle = angle
        self.unit = unit
        text_offset = 16

        self._used_deltas = ([], [])

        if minus_multi_marker is None:
            minus_multi_marker = MultiMarkerNode(
                mid_marker.x + self.unit * math.cos(self.angle + math.pi),
                mid_marker.y + self.unit * math.sin(self.angle + math.pi),
            )
        if minus_multi_marker is mid_marker:
            minus_multi_marker = None
        if plus_multi_marker is None:
            plus_multi_marker = MultiMarkerNode(
                mid_marker.x + self.unit * math.cos(self.angle),
                mid_marker.y + self.unit * math.sin(self.angle),
            )
        if plus_multi_marker is mid_marker:
            plus_multi_marker = None

        if label_x is None or label_y is None:
            label_x = mid_marker.x + text_offset * abs(math.cos(self.angle - 0.5 * math.pi))
            label_y = (
                mid_marker.y
                + text_offset * math.sin(self.angle - 0.5 * math.pi)
                + text_y_correction
            )

        super().__init__(
            bigg_id=bigg_id,
            mid_marker=mid_marker,
            label_x=label_x,
            label_y=label_y,
            minus_multi_marker=minus_multi_marker,
            plus_multi_marker=plus_multi_marker,
            **kwargs
        )

    def alternating_side_placement(self, i, delta, plus_minus):
        n = 1 + (i - 1) // 2
        side = (i - 1) % 2
        d = (n * delta) if side else -(n * delta)
        return n, side, d

    def same_side_placement(self, n, delta, plus_minus, absolute_side=0):
        side = bool(absolute_side) == bool(plus_minus)
        if (self.angle % (2 * math.pi)) > math.pi:
            side = not side
        side = int(side)

        d = (n * delta) if side else -(n * delta)
        return n, side, d

    def calculate_placement(
        self, ref_node, node, plus_minus, angle_delta, n, b1_b2, placement_opts
    ):
        size = placement_opts.scale
        angle = self.angle + angle_delta

        x, y = node.x, node.y
        if x is None or y is None:
            if not node.node_is_primary:
                size = placement_opts.no_primary_length_f(n) * placement_opts.scale
            x = ref_node.x + self.unit * size * math.cos(
                angle + (1 - plus_minus) * math.pi
            )
            y = ref_node.y + self.unit * size * math.sin(
                angle + (1 - plus_minus) * math.pi
            )
        else:
            size = math.sqrt((x - ref_node.x) ** 2 + (y - ref_node.y) ** 2) / self.unit

        if b1_b2 is not None:
            b1, b2 = b1_b2
        else:
            b2 = (
                x
                + self.unit
                * (1 - placement_opts.b2_scale)
                * size
                * math.cos(angle + (2 - plus_minus) * math.pi),
                y
                + self.unit
                * (1 - placement_opts.b2_scale)
                * size
                * math.sin(angle + (2 - plus_minus) * math.pi),
            )
            # b2 = None
            b1 = (
                ref_node.x
                + self.unit
                * placement_opts.b1_scale
                * size
                * math.cos(self.angle + (1 - plus_minus) * math.pi),
                ref_node.y
                + self.unit
                * placement_opts.b1_scale
                * size
                * math.sin(self.angle + (1 - plus_minus) * math.pi),
            )
        t = min(1.5 / size, 1.0)
        bt = (
            cubic_bezier_bt(t, ref_node.x, b1[0] if b1 is not None else ref_node.x, b2[0] if b2 is not None else x, x),
            cubic_bezier_bt(t, ref_node.y, b1[1] if b1 is not None else ref_node.y, b2[1] if b2 is not None else y, y),
        )
        effective_angle_delta = (
            math.atan2(bt[1] - ref_node.y, bt[0] - ref_node.x) - self.angle
        )
        if not plus_minus:
            effective_angle_delta = effective_angle_delta + math.pi
        effective_angle_delta = math.remainder(effective_angle_delta, math.pi * 2)

        return x, y, size, b1, b2, effective_angle_delta

    def add_metabolite(
        self,
        node: MetaboliteNode,
        coefficient: Union[float, int],
        b1_b2: Optional[Tuple[Optional[float], Optional[float]]] = None,
        placement_opts=None,
    ):
        if placement_opts is None:
            placement_opts = PlacementOptions()
        if placement_opts.placement_f is None:
            placement_opts.placement_f = self.__class__.alternating_side_placement
        ref_node = self.mid_marker
        plus_minus = coefficient > 0
        if self.multi_markers[plus_minus] is not None:
            ref_node = self.multi_markers[plus_minus]

        if node.x is not None and node.y is not None and b1_b2 is not None:
            x, y, size, b1, b2, effective_angle_delta = self.calculate_placement(
                ref_node, node, plus_minus, 0, 0, b1_b2, placement_opts
            )
            side = effective_angle_delta >= 0
            n = abs(effective_angle_delta) / placement_opts.delta
        else:
            for i in range(10):
                n, side, angle_delta = placement_opts.placement_f(
                    self, i, placement_opts.delta, plus_minus
                )
                x, y, size, b1, b2, effective_angle_delta = self.calculate_placement(
                    ref_node, node, plus_minus, angle_delta, n, b1_b2, placement_opts
                )
                # print(f"{node.bigg_id}: ({i}): (n:{n}, side:{side}, angle_delta:{angle_delta}) -> (x:{x}, y:{y}, size:{size}, b1:{b1}, b2:{b2}, effective_angle_delta:{effective_angle_delta})")
                if not any(
                    abs(d - effective_angle_delta) < (placement_opts.delta_tolerance * placement_opts.delta)
                    for d in self._used_deltas[plus_minus]
                ):
                    break

        self._used_deltas[plus_minus].append(effective_angle_delta)
        node.x = x
        node.y = y

        if node.label_x is None or node.label_y is None:
            x_positive = -min(
                0,
                math.cos(
                    self.angle
                    + (1 if bool(side) == bool(plus_minus) else -1) * 0.5 * math.pi
                ),
            )
            label_x = node.x + placement_opts.text_offset_f(
                x_positive * len(node.bigg_id)
            ) * math.cos(
                self.angle
                + (1 if bool(side) == bool(plus_minus) else -1) * 0.5 * math.pi
            )
            label_y = (
                node.y
                + placement_opts.text_offset_f(x_positive * len(node.bigg_id))
                * math.sin(
                    self.angle
                    + (1 if bool(side) == bool(plus_minus) else -1) * 0.5 * math.pi
                )
                + placement_opts.text_y_correction
            )
            node.label_x = label_x
            node.label_y = label_y

        self.metabolites.append((coefficient, node))
        if self._map is not None:
            self._map.add_node(node)

        if plus_minus:
            segment = Segment(ref_node, node, b1=b1, b2=b2)
        else:
            segment = Segment(node, ref_node, b1=b2, b2=b1)
        self.add_segment(segment)

class AutoReactionWithOptionalMetabolites(AutoReaction):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.optional_metabolites = {}
        self.finalized = False
    
    def add_optional_metabolite(self, node: MetaboliteNode, coefficient: float, b1_b2: Optional[Tuple[Optional[float], Optional[float]]]):
        self.optional_metabolites[node.bigg_id] = {"coefficient": coefficient, "node": node, "b1_b2": b1_b2}