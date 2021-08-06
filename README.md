# Nucleic-Transformer-WebApp



I made this web application (inspired by https://github.com/fatihozturkh2o/explorna_wave
) using h2o's wave so that models in my paper (https://www.biorxiv.org/content/10.1101/2021.01.28.428629v1) can be used easily. In addition, if your computer is too slow to run the app, you can also use Kaggle notebooks (you get 40 free GPU hours per week) I created to use the models:

https://www.kaggle.com/shujun717/nucleic-transformer-promoter-inference <br />
https://www.kaggle.com/shujun717/nucleic-transformer-virus-inference <br />
https://www.kaggle.com/shujun717/nucleic-transformer-rna-degradation-inference <br />

## Showcase

### Home page
![home_page](https://github.com/Shujun-He/Nucleic-Transformer-WebApp/blob/main/files/home_page.png)

### Promoter classification
Here you can classify DNA promoters and visualize the top kmers
![Promoter](https://github.com/Shujun-He/Nucleic-Transformer-WebApp/blob/main/files/promoter_page.png)

### Virus classification
Here you can classify DNA virus and visualize the top kmers
![Virus](https://github.com/Shujun-He/Nucleic-Transformer-WebApp/blob/main/files/virus_page.png)

### RNA degradation prediction
In this page you can predict RNA degradation at each nucleotide and visualize the attention weights of the Nucleic Transformer
![RNA degradation](https://github.com/Shujun-He/Nucleic-Transformer-WebApp/blob/main/files/rna_page.png)





## How to run
Clone the repo to your local. This should be run on a Linux Machine.

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
