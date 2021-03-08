git clone https://github.com/DasLab/arnie arnie;
git clone https://github.com/hendrixlab/bpRNA bpRNA;
cd arnie && git reset --hard bdb7f803fe24273c34b3402c3d18036d9bfc882e;
cd ..;
echo "TMP: tmp" > arnie.conf;
echo "rnasoft: MultiRNAFold" >> arnie.conf;
git clone https://www.github.com/DasLab/draw_rna draw_rna_pkg;
cd draw_rna_pkg && git reset --hard 90357a5ca529f834b50b07bc8f66ed198e43374d;
cd ..;
cp -R fixed_draw.py draw_rna_pkg/ipynb/draw.py;
cd draw_rna_pkg && python setup.py install;
cd ..
chmod u+x run.sh
cat Virus_Inference/best_weights/segmenta* > Virus_Inference/best_weights/fold0top1.ckpt
unzip MultiRNAFold.zip
