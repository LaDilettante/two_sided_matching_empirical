(TeX-add-style-hook
 "guide-notes"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("dukedissertation" "economy" "twoside" "bind")))
   (TeX-run-style-hooks
    "latex2e"
    "dukedissertation"
    "dukedissertation10")
   (LaTeX-add-labels
    "chap:guide"))
 :latex)

