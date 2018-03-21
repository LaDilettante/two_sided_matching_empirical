(TeX-add-style-hook
 "AnhLe_matching_progress_0304"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=1in")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "geometry"
    "graphicx"
    "pdfpages"
    "amsmath"
    "amsfonts"
    "amsthm"
    "bm"
    "enumitem"
    "rotating"
    "caption"
    "xcolor"
    "hyperref"
    "cleveref"
    "placeins"
    "float"
    "natbib")
   (LaTeX-add-bibliographies
    "/home/anh/Dropbox/texmf/bibtex/bib/library"))
 :latex)

