latexdiff v5352-20191029/main/octopus-paper.tex v5658-20200223/main/octopus-paper.tex  > v5658-20200223/main/octopus-paper-diff.tex
cd v5658-20200223/main
pdflatex octopus-paper-diff.tex
bibtex octopus-paper-diff
pdflatex octopus-paper-diff.tex
pdflatex octopus-paper-diff.tex
mv octopus-paper-diff.pdf ../..
rm octopus-paper-diff*
cd ../..

