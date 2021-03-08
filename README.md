# Nucleic-Transformer-WebApp

Based on https://github.com/fatihozturkh2o/explorna_wave

Clone the repo to your local.

Before everything make sure that wave is ready. 
[How to install]https://h2oai.github.io/wave/docs/installation/ Start the server with

./waved

This needs to keep running for the app to work. So leave the terminal with waved runinig and open a new one to do the rest of the setup

**1.** Download and install rnasoft from: http://www.rnasoft.ca/download/MultiRNAFold-2.1.tar.gz

To extract:
```bash
tar -xzvf MultiRNAFold-2.1.tar.gz
```

**2.** Open terminal in the cloned folder and run: <code>make setup</code>

**3.** Then run the bash file: for Ubuntu: <code>bash Ubuntu_setup_draw_rna.sh</code> , for Mac: <code>bash MacOS_setup_draw_rna.sh</code>

**4.** Now you are ready to run the app: <code>./venv/bin/python run.py</code>   

**5.** Go to <code>http://localhost:10101/nucleictransformer</code> in your browser.
