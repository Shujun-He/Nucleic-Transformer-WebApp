# Nucleic-Transformer-WebApp

Inspired by https://github.com/fatihozturkh2o/explorna_wave

Clone the repo to your local.

Before everything make sure that wave is ready. 
For instructions on installation: https://h2oai.github.io/wave/docs/installation/ 

Start the server with

```bash
cd $HOME/wave
./waved
```

This needs to be running for the app to work. So leave the terminal with waved running and open a new one to do the rest of the setup

**0.** Install Graph for perl. I use bpRNA to generate features so you will need perl to run it

```bash
yes '' | cpan -i Graph
```

**1.** Open terminal in the cloned folder and run: <code>make setup</code>

**2.** Then run the bash file: for Ubuntu: <code>bash Ubuntu_setup_draw_rna.sh</code> , for Mac: <code>bash MacOS_setup_draw_rna.sh</code>

**3.** Now you are ready to run the app: <code>./venv/bin/python run.py</code>   

**3.** Go to <code>http://localhost:10101/nucleictransformer</code> in your browser.
