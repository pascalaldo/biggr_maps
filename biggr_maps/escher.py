from escher import Builder
from IPython.display import IFrame, display
import tempfile

def display_builder(builder):
    tmp_filename = "/tmp/escher.html"
    builder.save_html(tmp_filename)

    iframe = IFrame(tmp_filename, width="100%", height="800")
    return iframe
    