from ipycanvas import Canvas

class IpycanvasDrawer:
    def __init__(self, *, width=None, height=None, save_image=None, **kwargs):
        self.canvas = Canvas(width=width, height=height, sync_image_data=True)

        if 'save_image' in kwargs and 0:
            # define callback per
            # https://ipycanvas.readthedocs.io/en/latest/retrieve_images.html
            def save_to_file(*args, **kwargs):
                self.canvas.to_file(self.save_image)

            self.canvas.observe(save_to_file, "image_data")

    def draw_rectangle(self, coords, color):
        self.canvas.fill_style = color
        self.canvas.fill_rect(*coords)

    def image(self):
        return self.canvas
