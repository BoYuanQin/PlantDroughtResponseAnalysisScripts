#!/bin/bash

muscle -super5 combined_sequences.fasta -output aligned_sequences_1.fasta -threads 128 

# 检查 MUSCLE 比对是否成功
if [ $? -ne 0 ]; then
    echo "MUSCLE 比对失败。请检查输入文件格式和 MUSCLE 安装。"
    exit 1
fi
#构建树 
nohup iqtree -s combined_sequences_extraction1_alignment.fasta -m LG+G+F -st AA -nt AUTO -bb 1000 -alrt 1000 -ninit 10 -ntop 5 -seed 12345 &

nohup iqtree -s combined_sequences_extraction2_alignment.fasta -m LG+G+F -st AA -nt AUTO -ninit 10 -ntop 5 -seed 12345 &