julia runPIDC.jl Input_Data/Expr.txt Output_Results/network_PIDC.txt

python runArboreto.py --algo GRNBoost2 --inFile Input_Data/Expr.txt --outFile Output_Results/Network_GRNBoost2.tsv

python runArboreto.py --algo GENIE3 --inFile Input_Data/Expr.txt --outFile Output_Results/Network_GENIE3.tsv