from PIL import Image, ImageDraw


class PillowDrawer:
    def __init__(self, *, width=None, height=None, **kwargs):
        self.im = Image.new('RGB', (width, height), 'white')
        self.imdraw = ImageDraw.Draw(self.im)

    def draw_rectangle(self, coords, color):
        x, y, width, height = coords
        self.imdraw.rectangle((x, y, x + width, y + height), fill=color)

    def image(self):
        return self.im
