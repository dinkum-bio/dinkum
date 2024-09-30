__all__ = ["MultiTissuePanel",
           "TissueActivityPanel" ]

from .draw_ipycanvas import IpycanvasDrawer
from .draw_pillow import PillowDrawer

class MultiTissuePanel:
    """
    Draw multiple panels, representing multiple tissues & gene expression
    in each one.
    """
    def __init__(self, *, states=None, tissue_names=None, save_image=None,
                 canvas_type=None, genes_by_name=None):
        # create individual panels for each tissue
        self.panels = [ TissueActivityPanel(states=states, tissue_name=t,
                                            genes_by_name=genes_by_name)
                        for t in tissue_names ]

        self.save_image = save_image
        self.canvas_type = canvas_type

    def draw(self, is_active_fn):
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
            x_offsets.append(width)

        # build the appropriate canvas type
        if self.canvas_type == 'ipycanvas':
            canvas = IpycanvasDrawer(width=total_width,
                                     height=max_height,
                                     save_image=self.save_image)
        else:
            canvas = PillowDrawer(width=total_width,
                                  height=max_height)

        # draw each tissue, with each collection of genes,
        # spread out horizontally
        for p, x_offset in zip(self.panels, x_offsets):
            # draw background
            d = p.draw_tissue(canvas, x_offset=x_offset)

            # draw time point/tissue/state
            gene_names = p.gene_names
            times = p.times
            d.draw(canvas, times, gene_names, is_active_fn)
            
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
    
    def __init__(self, *, states=None, tissue_name=None, genes_by_name=None):
        assert tissue_name is not None
        self.tissue_name = tissue_name

        assert states is not None

        # determine all genes relevant to this tissue, + times, from 'states'.
        times = []
        all_gene_names = set()
        for (tp, state) in states:
            times.append(tp)
            
            activity = state.get_by_tissue_name(tissue_name)
            all_gene_names.update(activity.genes_by_name)

        ordered_names = []
        if not genes_by_name:
            ordered_names = list(all_gene_names)
        else:
            # prioritize genes_by_name; do remainder alphabetically
            for k in genes_by_name:
                ordered_names.append(k)
            for k in sorted(all_gene_names):
                if k not in genes_by_name:
                    ordered_names.append(k)
        self.gene_names = ordered_names

        self.times = times
        self.states = states

    def estimate_panel_size(self):
        "Estimate the size of this panel, based on # times / # genes"
        height = (len(self.times) + 1) * (self.box_size + self.box_spacing) + \
            self.box_y_start
        width = len(self.gene_names) * (self.box_size + self.box_spacing) + \
            self.box_x_start
        return width, height
    
    def draw_tissue(self, canvas, *, x_offset=0):
        """Draw this tissue on existing canvas.

        Returns a TissueActivityPanel_Draw that can be used to fill in the
        actual gene/time point/tissue activity.
        """
        gene_names = self.gene_names
        times = self.times

        box_total_size = self.box_size + self.box_spacing

        locations_by_tg = {}

        # determine the locations of the time points / genes.
        for row in range(0, len(times)):
            for col in range(0, len(gene_names)):
                xpos = self.box_x_start + box_total_size*col + x_offset
                ypos = self.box_y_start + box_total_size*row
                
                coords = [(xpos, ypos),
                          (xpos + self.box_size, ypos),
                          (xpos + self.box_size, ypos + self.box_size),
                          (xpos, ypos + self.box_size)]

                timep = times[row]
                gene_name = gene_names[col]
                loc = Tissue_TimePointGene_Location(timep, gene_name, coords)

                # save!
                locations_by_tg[(timep, gene_name)] = loc

        tissue_label_ypos = len(gene_names) * box_total_size + \
            self.box_y_start
        tissue_label_xpos = self.box_x_start + x_offset + \
            round((box_total_size * len(gene_names)) / 2.0)

        canvas.draw_text(self.tissue_name, tissue_label_xpos,
                         tissue_label_ypos, align="center")
        self.locations_by_tg = locations_by_tg

        # draw row names / time points
        for row in range(0, len(times)):
            xpos = self.box_x_start - box_total_size / 2 + x_offset
            ypos = self.box_y_start + box_total_size*row
            canvas.draw_text(times[row], xpos, ypos, align="right")

        # draw col names / genes
        for col in range(0, len(gene_names)):
            ypos = self.box_y_start - box_total_size

            xpos = self.box_x_start + box_total_size*col
            xpos += x_offset

            canvas.draw_text(gene_names[col], xpos, ypos,
                             align="center")
                             #max_width = box_total_size)

        return TissueActivityPanel_Draw(self)


class TissueActivityPanel_Draw:
    "Use the timep/gene location to draw gene activity."
    active_color = "DeepSkyBlue"
    inactive_color = "DarkGrey"
    
    def __init__(self, template):
        self.template = template

    def draw(self, canvas, times, gene_names, is_active_fn):
        locations_by_tg = self.template.locations_by_tg
        tissue_name = self.template.tissue_name

        for tp in times:
            for gene in gene_names:
                color = None
                if is_active_fn(tissue_name, tp, gene):
                    color = self.active_color
                else:
                    color = self.inactive_color

                loc = locations_by_tg.get((tp, gene))
                if loc:
                    loc.draw(canvas, color)


###

class SeaUrchin_Blastula_ActivityPanel:
    def __init__(self, *, states=None, tissue_name=None):
        self.tissue_name = name
        self.states = states

    def draw_tissue(self, canvas, *, x_offset=0):
        pass
