__all__ = ["MultiTissuePanel",
           "TissueActivityPanel" ]

from ipycanvas import Canvas

class MultiTissuePanel:
    def __init__(self, *, states=None, tissue_names=None, save_image=None):
        self.panels = [ TissueActivityPanel(states=states, tissue_name=t) for t in tissue_names ]

        self.save_image = save_image
        
    def draw(self, is_active_fn):
        total_width = 0
        x_offsets = [0]
        max_height = 0
        for p in self.panels:
            width, height = p.estimate_panel_size()
            max_height = max(height, max_height)
            total_width += width
            x_offsets.append(width)

        canvas = Canvas(width=total_width, height=max_height,
                        sync_image_data=True)

        if self.save_image:
            # define callback per
            # https://ipycanvas.readthedocs.io/en/latest/retrieve_images.html
            def save_to_file(*args, **kwargs):
                canvas.to_file(self.save_image)

            canvas.observe(save_to_file, "image_data")
        
        for p, x_offset in zip(self.panels, x_offsets):
            p.draw_tissue(canvas, is_active_fn, x_offset=x_offset)
            
        return canvas

        
class TissueActivityPanel:
    box_size = 25
    box_spacing = 5
    
    box_x_start = 100
    box_y_start = 50
    
    active_color = "DeepSkyBlue"
    inactive_color = "DarkGrey"
    
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
    
    def draw(self, is_active):
        "Create canvas & draw single tissue."
        width, height = self.estimate_panel_size()
        canvas = Canvas(width=width, height=height,
                        sync_image_data=True)
        self.draw_tissue(canvas, is_active)
        
        return canvas

    def draw_tissue(self, canvas, is_active_fn, *, x_offset=0):
        "Draw this tissue on existing canvas."
        gene_names = self.gene_names
        times = self.times

        box_total_size = self.box_size + self.box_spacing

        for row in range(0, len(times)):
            for col in range(0, len(gene_names)):
                if is_active_fn(self.tissue_name, times[row], gene_names[col]):
                    canvas.fill_style = self.active_color
                else:
                    canvas.fill_style = self.inactive_color

                xpos = self.box_x_start + box_total_size*col + x_offset
                ypos = self.box_y_start + box_total_size*row
                canvas.fill_rect(xpos, ypos, self.box_size, self.box_size)

        canvas.font = "18px Arial"
        canvas.text_baseline = "top"
        canvas.fill_style = "black"

        # row names / time points
        canvas.text_align = "right"
        for row in range(0, len(times)):
            xpos = self.box_x_start - box_total_size / 2 + x_offset
            ypos = self.box_y_start + box_total_size*row
            canvas.fill_text(times[row], xpos, ypos)

        # col names / genes
        canvas.text_align = "center"
        for col in range(0, len(gene_names)):
            ypos = self.box_y_start - box_total_size
            xpos = self.box_x_start + box_total_size*col + box_total_size / 2 + x_offset

            canvas.fill_text(gene_names[col], xpos, ypos, max_width = box_total_size)