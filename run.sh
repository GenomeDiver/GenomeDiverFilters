#!/bin/bash 
java -Xms2g -Xmx4g --add-opens=java.base/sun.nio.ch=ALL-UNNAMED -jar ./target/scala-2.13/recommendation-assembly-0.1.0-SNAPSHOT.jar -e 20 -h HP:0001156,HP:0001363,HP:0010055,HP:0011304 -g /home/kshi/gd-data/patient/patient-1/analysis-1/rare.genes.tsv -l /home/kshi/gd-data/lirical-cli-2.0.0-RC2-json-old/data -o /home/kshi/gd-data/patient/patient-1/analysis-1/output.txt
