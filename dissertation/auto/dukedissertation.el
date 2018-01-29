(TeX-add-style-hook
 "dukedissertation"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("report" "12pt" "letterpaper")))
   (TeX-run-style-hooks
    "latex2e"
    "graphicx"
    "ifthen"
    "report"
    "rep12")
   (TeX-add-symbols
    '("dedication" 1)
    '("member" 1)
    '("copyrighttext" 1)
    '("supervisorb" 1)
    '("supervisora" 1)
    '("department" 1)
    "biographyname"
    "ackname"
    "loaname"
    "biography"
    "introduction"
    "abbreviations"
    "acknowledgements"
    "maketitle"
    "title"
    "author"
    "date"
    "and"
    "nmchapter")
   (LaTeX-add-environments
    "symbollist"))
 :latex)

