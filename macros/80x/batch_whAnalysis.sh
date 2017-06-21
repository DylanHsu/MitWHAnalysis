mkdir -p MitWHAnalysis/macros/80x
mv *.C MitWHAnalysis/macros/80x
mv *.h MitWHAnalysis/macros/80x
ls MitWHAnalysis/macros/80x
root -b -l -q MitWHAnalysis/macros/80x/whAnalysis.C+\(\"\",\"$3\",true,\"$1\",$2,$4,$5\)
