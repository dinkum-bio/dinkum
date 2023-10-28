__all__ = ["MultiTissuePanel",
           "TissueActivityPanel" ]

from .draw_ipycanvas import IpycanvasDrawer
from .draw_pillow import PillowDrawer

class MultiTissuePanel:
    def __init__(self, *, states=None, tissue_names=None, save_image=None,
                 canvas_type=None):
        self.panels = [ TissueActivityPanel(states=states, tissue_name=t) for t in tissue_names ]

        self.save_image = save_image
        self.canvas_type = canvas_type
        
    def draw(self, is_active_fn):
        total_width = 0
        x_offsets = [0]
        max_height = 0
        for p in self.panels:
            width, height = p.estimate_panel_size()
            max_height = max(height, max_height)
            total_width += width
            x_offsets.append(width)

        if self.canvas_type == 'ipycanvas':
            canvas = IpycanvasDrawer(width=total_width,
                                     height=max_height,
                                     save_image=self.save_image)
        else:
            canvas = PillowDrawer(width=total_width,
                                  height=max_height)

        for p, x_offset in zip(self.panels, x_offsets):
            # draw background
            d = p.draw_tissue(canvas, x_offset=x_offset)

            # draw time point/tissue/state
            gene_names = p.gene_names
            times = p.times
            d.draw(canvas, times, gene_names, is_active_fn)
            
        return canvas.image()


class Tissue_TimePointGene_Location:
    def __init__(self, timepoint, gene, polygon_coords):
        self.tp = timepoint
        self.gene = gene
        self.polygon_coords = polygon_coords

    def draw(self, canvas, color):
        canvas.draw_rectangle(self.polygon_coords, color)


class TissueActivityPanel:
    box_size = 25
    box_spacing = 5
    
    box_x_start = 100
    box_y_start = 50
    
    def __init__(self, *, states=None, tissue_name=None):
        assert tissue_name is not None
        self.tissue_name = tissue_name
    
        times = []
        all_gene_names = set()
        for (tp, state) in states:
            times.append(tp)
            
            activity = state.get_by_tissue_name(tissue_name)
            all_gene_names.update(activity.genes_by_name)
        self.gene_names = list(sorted(all_gene_names))

        self.times = times
        self.states = states

    def estimate_panel_size(self):
        height = len(self.times) * (self.box_size + self.box_spacing) + self.box_y_start
        width = len(self.gene_names) * (self.box_size + self.box_spacing) + self.box_x_start
        return width, height
    
    def draw_tissue(self, canvas, *, x_offset=0):
        "Draw this tissue on existing canvas."
        gene_names = self.gene_names
        times = self.times

        box_total_size = self.box_size + self.box_spacing

        locations_by_tg = {}

        for row in range(0, len(times)):
            for col in range(0, len(gene_names)):
                xpos = self.box_x_start + box_total_size*col + x_offset
                ypos = self.box_y_start + box_total_size*row
                
                coords = (xpos, ypos, self.box_size, self.box_size)
                loc = Tissue_TimePointGene_Location(times[row], gene_names[col], coords)
                locations_by_tg[(times[row], gene_names[col])] = loc

        self.locations_by_tg = locations_by_tg

        # row names / time points
        for row in range(0, len(times)):
            xpos = self.box_x_start - box_total_size / 2 + x_offset
            ypos = self.box_y_start + box_total_size*row
            canvas.draw_text(times[row], xpos, ypos, align="right")

        # col names / genes
        for col in range(0, len(gene_names)):
            ypos = self.box_y_start - box_total_size

            xpos = self.box_x_start + box_total_size*col
            xpos += x_offset

            canvas.draw_text(gene_names[col], xpos, ypos,
                             align="center")
                             #max_width = box_total_size)

        return TissueActivityPanel_Draw(self)


class TissueActivityPanel_Draw:
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
