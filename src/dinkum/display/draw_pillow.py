from PIL import Image, ImageDraw


class PillowDrawer:
    def __init__(self, *, width=None, height=None, **kwargs):
        self.im = Image.new('RGB', (width, height), 'white')
        self.imdraw = ImageDraw.Draw(self.im)

    def polygon(self, coords, fill):
        self.imdraw.polygon(coords, fill=fill)

    def draw_text(self, text, xpos, ypos, color="black", align="center"):
        if align == "center":
            anchor = "lm"
        elif align == "right":
            anchor = "rm"
        else:
            raise Exception(f"unknown align string: '{align}'")
            
        self.imdraw.text((xpos, ypos), text, fill=color, anchor=anchor)
        
    def image(self):
        return self.im
