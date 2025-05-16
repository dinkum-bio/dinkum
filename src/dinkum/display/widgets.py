__all__ = ["MultiTissuePanel", "TissueActivityPanel"]

from .draw_ipycanvas import IpycanvasDrawer
from .draw_pillow import PillowDrawer
from dinkum import vfg, vfn
import matplotlib

EMIT_WARNINGS = False

gene_cmap = matplotlib.colormaps.get_cmap("Blues")
receptor_off_cmap = matplotlib.colormaps.get_cmap("Reds")
receptor_on_cmap = matplotlib.colormaps.get_cmap("YlGn")
ligand_cmap = matplotlib.colormaps.get_cmap("Purples")


def map_to_color(level, is_receptor, is_ligand, is_active):
    if level < 0:
        print(f"warning: level is too small: {level}")
        level = 0
    if level > 100:
        print(f"warning: level is too big: {level}")
        level = 100

    f = (level / 100) * 0.6 + 0.2  # pick out the middle 60%
    if is_ligand:
        color = ligand_cmap(f)
    elif is_receptor:
        if is_active:
            color = receptor_on_cmap(f)
        else:
            color = receptor_off_cmap(f)
    else:
        color = gene_cmap(f)
    color = tuple([int(x * 255) for x in color])
    return color


class MultiTissuePanel:
    """
    Draw multiple panels, representing multiple tissues & gene expression
    in each one.
    """

    def __init__(
        self,
        *,
        states=None,
        tissue_names=None,
        save_image=None,
        canvas_type="pillow",
        gene_names=None,
    ):
        # create individual panels for each tissue
        self.panels = [
            TissueActivityPanel(states=states, tissue_name=t, gene_names=gene_names)
            for t in tissue_names
        ]

        self.save_image = save_image
        self.canvas_type = canvas_type

    def draw(self, gene_states, level=100):
        """
        Draw the basic background canvas, upon which gene activation info
        will be displayed.
        """

        # determine overall panel size
        total_width = 0
        x_offsets = [0]
        max_height = 0

        for p in self.panels:
            width, height = p.estimate_panel_size()
            max_height = max(height, max_height)
            total_width += width
            x_offsets.append(total_width)

        # build the appropriate canvas type
        if self.canvas_type == "ipycanvas":
            canvas = IpycanvasDrawer(
                width=total_width, height=max_height, save_image=self.save_image
            )
        else:
            assert self.canvas_type == "pillow"
            canvas = PillowDrawer(
                width=total_width, height=max_height, save_image=self.save_image
            )

        # draw each tissue, with each collection of genes,
        # spread out horizontally
        for p, x_offset in zip(self.panels, x_offsets):
            # draw background
            d = p.draw_tissue(canvas, x_offset=x_offset)

            # draw time point/tissue/state
            gene_names = p.gene_names
            d.draw(canvas, gene_states.keys(), gene_names, gene_states)

        if self.save_image:
            canvas.save()

        return canvas.image()


class Tissue_TimePointGene_Location:
    """
    A class to track each time point & gene location within a particular
    tissue.
    """

    def __init__(self, timepoint, gene, polygon_coords):
        self.tp = timepoint
        self.gene = gene
        self.polygon_coords = polygon_coords

    def draw(self, canvas, color):
        canvas.polygon(self.polygon_coords, fill=color)


class TissueActivityPanel:
    """
    A class to draw and display gene activity in a single tissue.

    Constructed to be used within a MultiTissuePanel, among other things.
    """

    box_size = 25
    box_spacing = 5

    box_x_start = 100
    box_y_start = 50

    def __init__(self, *, states=None, tissue_name=None, gene_names=None):
        assert tissue_name is not None
        self.tissue_name = tissue_name

        assert states is not None

        # determine all genes relevant to this tissue, + times, from 'states'.
        times = []
        times_str = []
        all_gene_names = set()
        for timepoint, state in states.items():
            tp_str = f"t={timepoint}"
            times.append(timepoint)
            times_str.append(tp_str)

            activity = state.get_by_tissue_name(tissue_name)
            all_gene_names.update(activity.genes_by_name)

        ordered_names = []
        if not gene_names:
            ordered_names = list(all_gene_names)
        else:
            # prioritize gene_names; do remainder alphabetically
            for k in gene_names:
                ordered_names.append(k)
        #            for k in sorted(all_gene_names):
        #                if k not in gene_names:
        #                    ordered_names.append(k)
        self.gene_names = ordered_names

        self.times = times
        self.times_str = times_str
        self.states = states

    def estimate_panel_size(self):
        "Estimate the size of this panel, based on # times / # genes"
        height = (len(self.times_str) + 1) * (
            self.box_size + self.box_spacing
        ) + self.box_y_start
        width = (
            len(self.gene_names) * (self.box_size + self.box_spacing) + self.box_x_start
        )
        return width, height

    def draw_tissue(self, canvas, *, x_offset=0):
        """Draw this tissue on existing canvas.

        Returns a TissueActivityPanel_Draw that can be used to fill in the
        actual gene/time point/tissue activity.
        """
        gene_names = self.gene_names
        times = self.times
        times_str = self.times_str

        box_total_size = self.box_size + self.box_spacing

        locations_by_tg = {}

        # determine the locations of the time points / genes.
        for row in range(0, len(times)):
            for col in range(0, len(gene_names)):
                xpos = self.box_x_start + box_total_size * col + x_offset
                ypos = self.box_y_start + box_total_size * row

                coords = [
                    (xpos, ypos),
                    (xpos + self.box_size, ypos),
                    (xpos + self.box_size, ypos + self.box_size),
                    (xpos, ypos + self.box_size),
                ]

                timep = times[row]
                gene_name = gene_names[col]
                loc = Tissue_TimePointGene_Location(timep, gene_name, coords)

                # save!
                locations_by_tg[(timep, gene_name)] = loc

        tissue_label_ypos = (
            len(times) * box_total_size + self.box_y_start + 1 / 2 * box_total_size
        )
        tissue_label_xpos = (
            self.box_x_start
            + x_offset
            + round((box_total_size * len(gene_names)) / 2.0)
        )

        canvas.draw_text(
            self.tissue_name, tissue_label_xpos, tissue_label_ypos, align="center"
        )
        self.locations_by_tg = locations_by_tg

        # draw row names / time points
        for row in range(0, len(times_str)):
            xpos = self.box_x_start - box_total_size / 2 + x_offset
            ypos = self.box_y_start + box_total_size * row
            canvas.draw_text(times_str[row], xpos, ypos, align="right")

        # draw col names / genes
        for col in range(0, len(gene_names)):
            ypos = self.box_y_start - box_total_size

            xpos = self.box_x_start + box_total_size * col
            xpos += x_offset

            canvas.draw_text(gene_names[col], xpos, ypos, align="center")
            # max_width = box_total_size)

        return TissueActivityPanel_Draw(self)


class TissueActivityPanel_Draw:
    "Use the timep/gene location to draw gene activity."

    active_color = "DeepSkyBlue"  # blue!
    active_receptor_color = "DarkOliveGreen"  # blue!

    present_color = (255, 0, 0)  # red!
    present_mask = (0, 255, 255)  # ...?
    inactive_color = "DarkGrey"  # grey...

    def __init__(self, template):
        self.template = template

    def draw(self, canvas, times, gene_names, gene_states):
        locations_by_tg = self.template.locations_by_tg
        tissue_name = self.template.tissue_name
        tissue_obj = vfn.get_tissue(tissue_name)

        for tp in times:
            for gene_name in gene_names:
                color = None

                active_color = self.active_color
                gene_obj = vfg.get_gene(gene_name)
                gsi = gene_states.get_gene_state_info(
                    timepoint=tp, delay=0, gene=gene_obj, tissue=tissue_obj
                )

                # check a few constraints...
                if gsi.level == 0 and gsi.active and EMIT_WARNINGS:
                    print(
                        f"WARNING: at time {tp} gene {gene_name} has level == 0 but is active! Is this intentional??"
                    )
                if (
                    not gene_obj.is_receptor
                    and gsi.level > 0
                    and not gsi.active
                    and EMIT_WARNINGS
                ):
                    print(
                        f"WARNING: at time {tp} gene {gene_name} has level {gsi.level} but is not active! Is this intentional??"
                    )
                    raise Exception("what 2")

                color = map_to_color(
                    gsi.level, gene_obj.is_receptor, gene_obj.is_ligand, gsi.active
                )

                loc = locations_by_tg.get((tp, gene_name))
                if loc:
                    loc.draw(canvas, color)


###


class SeaUrchin_Blastula_ActivityPanel:
    def __init__(self, *, states=None, tissue_name=None):
        self.tissue_name = name
        self.states = states

    def draw_tissue(self, canvas, *, x_offset=0):
        pass
