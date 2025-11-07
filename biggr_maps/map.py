import math
from typing import Dict, List, Optional, Tuple, Union


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
        self.text_labels = {}
        self.canvas = canvas
        self._node_counter = 0
        self._reaction_counter = 0
        self._segment_counter = 0

    def add_node(self, node: Optional["Node"]):
        if node is None or node.identifier is not None:
            return
        node.identifier = self._node_counter
        self.nodes[node.identifier] = node
        self._node_counter += 1

    def add_reaction(self, reaction: "Reaction"):
        reaction.identifier = self._reaction_counter
        self.reactions[reaction.identifier] = reaction
        self._reaction_counter += 1

        for _, node in reaction.metabolites:
            self.add_node(node)
        self.add_node(reaction.mid_marker)
        self.add_node(reaction.multi_markers[0])
        self.add_node(reaction.multi_markers[1])

        for segment in reaction.segments:
            segment.identifier = self._segment_counter
            self._segment_counter += 1

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

    def _add_multi_marker_segments(self):
        if self.multi_markers[0] is not None:
            self.segments.append(Segment(self.multi_markers[0], self.mid_marker))
        if self.multi_markers[1] is not None:
            self.segments.append(Segment(self.mid_marker, self.multi_markers[1]))

    def add_metabolite(self, node: MetaboliteNode, coefficient: Union[float, int]):
        self.metabolites.append((coefficient, node))

        ref_node = self.mid_marker
        plus_minus = coefficient > 0
        if self.multi_markers[plus_minus] is not None:
            ref_node = self.multi_markers[plus_minus]

        if plus_minus:
            segment = Segment(ref_node, node)
        else:
            segment = Segment(node, ref_node)
        self.segments.append(segment)

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


class AutoReaction(Reaction):
    def __init__(
        self,
        mid_marker: MidMarkerNode,
        angle: float,
        unit: float=50,
        text_y_correction: float=8,
        **kwargs
    ):
        self.angle = angle
        self.unit = unit
        text_offset = 16

        self._used_deltas = ([], [])

        multi1 = MultiMarkerNode(
            mid_marker.x + self.unit * math.cos(self.angle + math.pi),
            mid_marker.y + self.unit * math.sin(self.angle + math.pi),
        )
        multi2 = MultiMarkerNode(
            mid_marker.x + self.unit * math.cos(self.angle),
            mid_marker.y + self.unit * math.sin(self.angle),
        )

        label_x = mid_marker.x + abs(text_offset * math.cos(self.angle - 0.5 * math.pi))
        label_y = (
            mid_marker.y
            + text_offset * math.sin(self.angle - 0.5 * math.pi)
            + text_y_correction
        )

        super().__init__(
            mid_marker=mid_marker,
            label_x=label_x,
            label_y=label_y,
            minus_multi_marker=multi1,
            plus_multi_marker=multi2,
            **kwargs
        )

    def _determine_n_and_side(self, desired_delta, delta, plus_minus):
        if desired_delta is not None:
            if not any(
                abs(x - desired_delta) < delta for x in self._used_deltas[plus_minus]
            ):
                side = desired_delta >= 0
                n = abs(desired_delta) / delta
                return n, side, desired_delta
        i = 1
        while True:
            n = 1 + (i - 1) // 2
            side = (i - 1) % 2
            d = (n * delta) if side else -(n * delta)
            if not any(abs(x - d) < delta for x in self._used_deltas[plus_minus]):
                return n, side, d
            i += 1

    def add_metabolite(
        self,
        node,
        coefficient: Union[float, int],
        scale: float = 3.0,
        delta=math.pi * 0.15,
        no_primary_length_f=lambda x: (1 - (min(x - 1, 5) / 5)) * 0.3 + 0.5,
        b1_scale=0.3,
        b2_scale=0.8,
        text_y_correction=6,
        text_offset_f=lambda x: 20 + x * 12,
        b1_b2: Optional[Tuple[Optional[float], Optional[float]]] = None,
    ):
        ref_node = self.mid_marker
        plus_minus = coefficient > 0
        if self.multi_markers[plus_minus] is not None:
            ref_node = self.multi_markers[plus_minus]

        desired_delta = None
        if node.x is not None and node.y is not None:
            dy = node.y - ref_node.y
            dx = node.x - ref_node.x
            desired_delta = math.atan2(dy, dx) - self.angle - (1 - plus_minus) * math.pi

        is_primary = node.node_is_primary
        if desired_delta is None and is_primary:
            desired_delta = 0

        size = scale
        n, side, angle_delta = self._determine_n_and_side(
            desired_delta, delta, plus_minus
        )
        angle = self.angle + angle_delta
        self._used_deltas[plus_minus].append(angle_delta)
        bend_curve = abs(angle_delta) > 0.001
        if node.x is None or node.y is None:
            if not is_primary:
                size = no_primary_length_f(n) * scale
            x = ref_node.x + self.unit * size * math.cos(
                angle + (1 - plus_minus) * math.pi
            )
            y = ref_node.y + self.unit * size * math.sin(
                angle + (1 - plus_minus) * math.pi
            )
            node.x = x
            node.y = y
        else:
            size = (
                math.sqrt((node.x - ref_node.x) ** 2 + (node.y - ref_node.y) ** 2)
                / self.unit
            )

        if node.label_x is None or node.label_y is None:
            x_positive = -min(
                0,
                math.cos(
                    self.angle
                    + (1 if bool(side) == bool(plus_minus) else -1) * 0.5 * math.pi
                ),
            )
            label_x = node.x + text_offset_f(x_positive * len(node.bigg_id)) * math.cos(
                self.angle
                + (1 if bool(side) == bool(plus_minus) else -1) * 0.5 * math.pi
            )
            label_y = (
                node.y
                + text_offset_f(x_positive * len(node.bigg_id))
                * math.sin(
                    self.angle
                    + (1 if bool(side) == bool(plus_minus) else -1) * 0.5 * math.pi
                )
                + text_y_correction
            )
            node.label_x = label_x
            node.label_y = label_y

        self.metabolites.append((coefficient, node))

        if b1_b2 is not None:
            b1, b2 = b1_b2
        elif bend_curve:
            b2 = (
                node.x
                + self.unit
                * (1 - b2_scale)
                * size
                * math.cos(angle + (2 - plus_minus) * math.pi),
                node.y
                + self.unit
                * (1 - b2_scale)
                * size
                * math.sin(angle + (2 - plus_minus) * math.pi),
            )
            b1 = (
                ref_node.x
                + self.unit
                * b1_scale
                * size
                * math.cos(self.angle + (1 - plus_minus) * math.pi),
                ref_node.y
                + self.unit
                * b1_scale
                * size
                * math.sin(self.angle + (1 - plus_minus) * math.pi),
            )

        if plus_minus:
            segment = Segment(ref_node, node, b1=b1, b2=b2)
        else:
            segment = Segment(node, ref_node, b1=b2, b2=b1)
        self.segments.append(segment)
