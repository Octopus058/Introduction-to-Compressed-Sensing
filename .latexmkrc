# fontspec needs a Unicode engine. We use LuaLaTeX.
# Setting this makes a bare `latexmk` (no flags) do the right thing.
$pdf_mode = 4;    # 4 = lualatex

# Without this, a bare `latexmk` builds every .tex in the directory, which is
# not what we want here. Restrict it to the book.
@default_files = ("Introduction to Compressed Sensing.tex");
