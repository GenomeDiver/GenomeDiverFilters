#!/bin/bash 
java -Xms2g -Xmx4g --add-opens=java.base/sun.nio.ch=ALL-UNNAMED -jar ./target/scala-2.13/recommendation-assembly-0.1.0-SNAPSHOT.jar
