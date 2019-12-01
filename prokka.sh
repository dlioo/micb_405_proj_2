#!/bin/bash

cd ~
if [ ! -d project_2 ]; then
  mkdir project_2
fi
cd ~/project_2
if [ ! -d prokka ]; then
  mkdir prokka
fi

mags_path="/projects/micb405/resources/project_2/2019/SaanichInlet_100m/MetaBAT2_SaanichInlet_100m/MedQPlus_MAGs"

for fa in ${mags_path}/*fa
do
  bin="$(echo $fa | cut -d'.' -f2 )"
   prokka \
   --outdir ~/project_2/prokka/$bin \
   --prefix "${bin}_SaanichInlet_100m_bacteria" \
   --kingdom Bacteria \
   --cpus 8 \
   --force \
   $fa

   bin="$(echo $fa | cut -d'.' -f2)"
    prokka \
    --outdir ~/project_2/prokka/$bin \
    --prefix "${bin}_SaanichInlet_100m_archaea" \
    --kingdom Archaea \
    --cpus 8 \
    --force \
    $fa
done
