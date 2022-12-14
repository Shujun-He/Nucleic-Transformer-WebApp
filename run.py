# Form / Frame
# Use a frame component in a form card to display HTML content inline.
# ---
import sys
sys.path.append('draw_rna_pkg/')
from h2o_wave import Q, listen, ui
import matplotlib.pyplot as plt
import os
import io
import base64
import matplotlib.image as mpimg
import zipfile
import time



os.environ["ARNIEFILE"] = f"arnie.conf"


from dna_analysis import *
import Promoter_Inference
import Virus_Inference
import Enhancer_Inference

#from model_prediction import load_models,bpps_path,get_prediction_df_dict,bpps_check






# #for cmaps check out https://matplotlib.org/3.2.1/tutorials/colors/colormaps.html
target_columns = [ 'reactivity', 'deg_Mg_pH10',
       'deg_pH10', 'deg_Mg_50C', 'deg_50C']

logo_file = "files/logo.png"
image_1_path = "files/NT_architecture.png"
image_2_path = "files/image_2_path.png"
image_3_path = "files/image_3_path.png"
image_4_path = "files/image_4_path.png"

all_pages = ['nav','home', 'upload', 'example', 'infoboard1', 'image0', 'image1', 'image2','data_info','stn_hist',
                 'image3', 'image4', 'image5', 'description', 'text','text2', 'error', 'progress', 'selection_card', 'plot',
                 "plot_seqpos","plot_count",'plot_features0','plot_features1','plot_features2','plot_features3','text_sub','predict_tab','predict_promoter_tab']


image_size = 6
number_of_plots_in_a_row = 3
plot_height = 5
plot_width = 3

data_display_max_nrows = 10



def get_image(file_name = "image1"):
    plt.figure(figsize=(image_size, image_size))
    buf = io.BytesIO()

    img = mpimg.imread(f'temp/{file_name}.png')
    plt.imshow(img)
    plt.axis('off')
    plt.savefig(buf, format='png')
    buf.seek(0)
    image = base64.b64encode(buf.read()).decode('utf-8')

    return image

def display_error(error_type):
    items = []
    if error_type == "upload":
        items = [ui.message_bar(type='error', text='**Selected file was not in the expected format.**'),
                    ui.button(name='upload_data', label='Upload again', primary=True),
                    ui.separator()]
    elif error_type == "misentered_info":
        items = [ui.message_bar(type='error', text='**Either sequence or structure info not entered.**'),
                ui.button(name='custom_back', label='Try again', primary=True),
                ui.separator()]
    elif error_type == "wrong_columns":
        items = [
                ui.message_bar(type='error', text='**Either sequence or structure column not selected properly**'),
                ui.button(name='restart', label='Try again', primary=True),
                ui.separator()
            ]
    elif error_type == "misentered_ids":
        items = [
                ui.message_bar(type='error', text='**IDs not entered properly**'),
                ui.button(name='restart', label='Try again', primary=True),
                ui.separator()
            ]
    elif error_type == "arnie_not_selected":
        items = [
            ui.message_bar(type='error', text='**Either feature or package not selected properly.**'),
            ui.button(name='arnie_back', label='Try again', primary=True),
            ui.separator()
        ]
    return items

def make_ui_table(df, n_rows, name = "head_of_table"):
    """Creates a ui.table object from a csv file"""

    n_rows = min(n_rows, df.shape[0])

    table = ui.table(
        name=name,
        columns=[ui.table_column(name=str(x), label=str(x), sortable=True) for x in df.columns.values],
        rows=[ui.table_row(name=str(i), cells=[str(df[col].values[i]) for col in df.columns.values])
              for i in range(n_rows)]
    )
    return table

async def delete_pages(q: Q,keep_nav=False):
    # There does not seem to be a way to know what pages are in q.page
    page_list = ['nav','home', 'upload', 'example', 'infoboard1','data_info','stn_hist',
                 'description', 'text','text2', 'error', 'progress', 'selection_card', 'plot','data_view',
                 "plot_seqpos","plot_count",'plot_features0','plot_features1','plot_features2','plot_features3','text_sub','predict_tab','predict_promoter_tab'] \
                + list(set(q.client.all_pages))
    if keep_nav:
        page_list.remove("nav")

    for page in page_list:
        try:
            del q.page[page]
        except:
            #print('some error')
            pass

async def display_nav(q: Q):

    q.page['nav'] = ui.tab_card(
        box='3 1 9 1',
        items=[
            ui.tab(name='#home', label='Home'),
            ui.tab(name='#promoter_prediction', label='Promoter Classification'),
            ui.tab(name='#enhancer_prediction', label='Enhancer Classification'),
            ui.tab(name='#virusprediction', label='Virus Prediction'),
        ],
        link=False
    )


    q.page['header2'] = ui.header_card(
        box='1 1 2 1',
        title='Nucleic Transformer',
        subtitle='',
        icon='ExploreData',
        icon_color='$orange',
    )

    q.page['logo'] = ui.markup_card(
        box='12 1 1 1',
        title='',
        content="""<p style='text-align:center; vertical-align: middle; display: table-cell; width: 134px;'>"""
                """<a href='https://www.h2o.ai/products/h2o-wave/'> <img src='""" + q.app.logo_url + """' height='50px' width='50px'> </a> </p>""" )

async def progress_page(q: Q, message='\nJust a second...'):
    q.page['progress'] = ui.form_card(box='4 5 6 2',
                                      items=[ui.progress(label=message)])
    await q.page.save()
    del q.page['progress']



async def display_file_upload(q):
    print("location: display_file_upload")

    q.page['description'] = ui.form_card(box='1 2 3 10',
        items=[ui.text_m('\n You can upload a local dataset to run ExploRNA.'),
               ui.file_upload(name='user_files', label='Upload', multiple=False)])

    if q.client.train is not None:
        data_items = [ui.text_m(f'Loaded file "{q.client.file_name}" has '
                                f'**{q.client.train.shape[0]}** rows and **{q.client.train.shape[1]}** features.\n\n'),
                      make_ui_table(q.client.train, data_display_max_nrows)]
        q.page['data_view'] = ui.form_card(box='4 2 9 7', items=data_items)





async def promoter_model_predict(q):

    # Loading models
    start_time_main = time.time()

    if q.client.models_loaded is None:
        q.page['progress'] = ui.form_card(box='4 6 9 1',
                                          items=[ui.progress(label="Loading models...")])
        await q.page.save()
        del q.page['progress']
        q.client.promoter_inference=Promoter_Inference.Promoter_Inference()
        q.client.promoter_inference.load_models('Promoter_Inference/best_weights')
        q.client.models_loaded=True

    # Getting predictions
    #if q.client.predictions is None:
    minimum_required_features = ["sequence"]
    #await progress_page(q,message="Getting predictions...")

    q.page['progress'] = ui.form_card(box='4 6 9 1',
                                      items=[ui.progress(label="Getting predictions...")])
    await q.page.save()
    del q.page['progress']

    q.client.promoter_predictions, q.client.top_kmers, q.client.top_kmer_counts = \
    q.client.promoter_inference.predict(q.client.promoter_data.loc[:, minimum_required_features])
    q.client.promoter_predictions.to_csv('temp/promoter_predictions.csv', index=False)
    q.client.predictions_path, = await q.site.upload(['temp/promoter_predictions.csv'])

    elapsed_time = time.time() - start_time_main
    print(f"minutes passed for getting predictions: {round(elapsed_time/60, 2)}")




async def predict_promoter_tab(q):
    if q.client.promoter_topk is not None:
        k_value=q.client.promoter_topk
    else:
        k_value=10

    if q.client.promoter_data is None:
        q.client.promoter_data=q.client.train



    #if q.client.predictions is not None:
    q.page['predict_promoter_tab'] = ui.form_card(box='1 2 3 11', items=[
        ui.text_m(f'In this section, you can classify DNA promoters'),
        ui.text_m('You can upload a local dataset to classify with the Nucleic Transformer. You need to put \
        the sequences\
        into a csv file with \'sequence\' as a column. We have uploaded a sample file for you so you can follow \
        the format easily'),
        ui.file_upload(name='promoter_user_files', label='Upload', multiple=False),
        ui.text_m('Aside from classifying promoters, you can also visualize the top kmers extracted based on \
        Nucleic Trasnformer\'s self-attention weights'),
        ui.slider(name='promoter_topk', label='Select number of top kmers to visualize', min=3, max=10, step=1, value=k_value),
        ui.button(name='predict_promoter', label='Predict', primary=True)
    ])


    #if q.client.promoter_data is not None:

    data_items = [ui.text_m(f'Loaded file "{q.client.file_name}" has '
                            f'**{q.client.promoter_data.shape[0]}** rows and **{q.client.promoter_data.shape[1]}** features.\n\n'),
                  make_ui_table(q.client.promoter_data, data_display_max_nrows)]
    q.page['promoter_data_view'] = ui.form_card(box='4 2 9 4', items=data_items)

    if 'promoter_data_view' not in q.client.all_pages:
        q.client.all_pages.append('promoter_data_view')

    if q.client.promoter_predictions is not None:

        download_data_text = '''=
Inference complete! Click [here]({{predictions}}) to download the predictions! Look below for visualization of promoter composition and top kmers!
'''
        q.page['download_promoter_predictions'] = ui.markdown_card(
                box='4 6 9 1',
                title='',
                content=download_data_text,
                data=dict(predictions=q.client.predictions_path)
            )
        q.client.promoter_data_predictions_seq_length= plot_promoter_percent(q.client.promoter_predictions)

        q.page['plot_promoters'] = ui.frame_card(
            box='4 7 4 6',
            title='How many promoters/non-promoters are in the dataset',
            content=q.client.promoter_data_predictions_seq_length
        )

        if 'download_promoter_predictions' not in q.client.all_pages:
            q.client.all_pages.append('download_promoter_predictions')

        if 'plot_promoters' not in q.client.all_pages:
            q.client.all_pages.append('plot_promoters')

        q.client.promoter_topk_plot=plot_top_promoter_kmers(q.client.top_kmers, q.client.top_kmer_counts, q.client.promoter_topk)



        q.page['plot_promoter_kmers'] = ui.frame_card(
            box='8 7 5 6',
            title='Top promoter kmers extracted from attention weights',
            content=q.client.promoter_topk_plot
        )


        if 'plot_promoter_kmers' not in q.client.all_pages:
            q.client.all_pages.append('plot_promoter_kmers')


async def enhancer_model_predict(q):

    # Loading models
    start_time_main = time.time()

    if q.client.enhancer_models_loaded is None:
        q.page['progress'] = ui.form_card(box='4 6 9 1',
                                          items=[ui.progress(label="Loading models...")])
        await q.page.save()
        del q.page['progress']
        q.client.enhancer_inference=Enhancer_Inference.Enhancer_Inference()
        q.client.enhancer_inference.load_models('Enhancer_Inference/best_weights')
        q.client.enhancer_models_loaded=True

    # Getting predictions
    #if q.client.predictions is None:
    minimum_required_features = ["sequence"]
    #await progress_page(q,message="Getting predictions...")

    q.page['progress'] = ui.form_card(box='4 6 9 1',
                                      items=[ui.progress(label="Getting predictions...")])
    await q.page.save()
    del q.page['progress']

    q.client.enhancer_predictions, q.client.top_kmers, q.client.top_kmer_counts = \
    q.client.enhancer_inference.predict(q.client.enhancer_data.loc[:, minimum_required_features])
    q.client.enhancer_predictions.to_csv('temp/enhancer_predictions.csv', index=False)
    q.client.predictions_path, = await q.site.upload(['temp/enhancer_predictions.csv'])

    elapsed_time = time.time() - start_time_main
    print(f"minutes passed for getting predictions: {round(elapsed_time/60, 2)}")


async def predict_enhancer_tab(q):
    if q.client.enhancer_topk is not None:
        k_value=q.client.enhancer_topk
    else:
        k_value=10

    if q.client.enhancer_data is None:
        q.client.enhancer_file="enhancer_sample.csv"
        q.client.enhancer_data = pd.read_csv('Enhancer_Inference/'+q.client.enhancer_file)
        q.client.enhancer_data["sequence_length"] = q.client.enhancer_data["sequence"].apply(lambda seq: len(seq))



    #if q.client.predictions is not None:
    q.page['predict_enhancer_tab'] = ui.form_card(box='1 2 3 11', items=[
        ui.text_m(f'In this section, you can classify DNA enhancers'),
        ui.text_m('You can upload a local dataset to classify with the Nucleic Transformer. You need to put \
        the sequences\
        into a csv file with \'sequence\' as a column. We have uploaded a sample file for you so you can follow \
        the format easily'),
        ui.file_upload(name='enhancer_user_files', label='Upload', multiple=False),
        ui.text_m('Aside from classifying enhancers, you can also visualize the top kmers extracted based on \
        Nucleic Trasnformer\'s self-attention weights'),
        ui.slider(name='enhancer_topk', label='Select number of top kmers to visualize', min=3, max=10, step=1, value=k_value),
        ui.button(name='predict_enhancer', label='Predict', primary=True)
    ])


    #if q.client.enhancer_data is not None:

    data_items = [ui.text_m(f'Loaded file "{q.client.enhancer_file}" has '
                            f'**{q.client.enhancer_data.shape[0]}** rows and **{q.client.train.shape[1]}** features.\n\n'),
                  make_ui_table(q.client.enhancer_data, data_display_max_nrows)]
    q.page['enhancer_data_view'] = ui.form_card(box='4 2 9 4', items=data_items)

    if 'enhancer_data_view' not in q.client.all_pages:
        q.client.all_pages.append('enhancer_data_view')

    if q.client.enhancer_predictions is not None:

        download_data_text = '''=
Inference complete! Click [here]({{predictions}}) to download the predictions! Look below for visualization of enhancer composition and top kmers!
'''
        q.page['download_enhancer_predictions'] = ui.markdown_card(
                box='4 6 9 1',
                title='',
                content=download_data_text,
                data=dict(predictions=q.client.predictions_path)
            )
        q.client.enhancer_data_predictions_seq_length= plot_enhancer_percent(q.client.enhancer_predictions)

        q.page['plot_enhancers'] = ui.frame_card(
            box='4 7 4 6',
            title='How many enhancers/non-enhancers are in the dataset',
            content=q.client.enhancer_data_predictions_seq_length
        )

        if 'download_enhancer_predictions' not in q.client.all_pages:
            q.client.all_pages.append('download_enhancer_predictions')

        if 'plot_enhancers' not in q.client.all_pages:
            q.client.all_pages.append('plot_enhancers')

        q.client.enhancer_topk_plot=plot_top_enhancer_kmers(q.client.top_kmers, q.client.top_kmer_counts, q.client.enhancer_topk)



        q.page['plot_enhancer_kmers'] = ui.frame_card(
            box='8 7 5 6',
            title='Top enhancer kmers extracted from attention weights',
            content=q.client.enhancer_topk_plot
        )


        if 'plot_enhancer_kmers' not in q.client.all_pages:
            q.client.all_pages.append('plot_enhancer_kmers')


async def predict_virus_tab(q):
    if q.client.virus_topk is not None:
        k_value=q.client.virus_topk
    else:
        k_value=5

    if q.client.virus_data is None:
        q.client.virus_data=q.client.virus_data



    #if q.client.predictions is not None:
    q.page['predict_virus_tab'] = ui.form_card(box='1 2 3 11', items=[
        ui.text_m(f'In this section, you can classify DNA virus'),
        ui.text_m('You can upload a local dataset to classify with the Nucleic Transformer. You need to put \
        the sequences\
        into a csv file with \'sequence\' as a column. We have uploaded a sample file for you so you can follow \
        the format easily'),
        ui.file_upload(name='virus_user_files', label='Upload', multiple=False),
        ui.text_m('In addition, you can also visualize the top kmers extracted based on Nucleic Trasnformer\'s self-attention weights'),
        ui.slider(name='virus_topk', label='Select number of top kmers to visualize', min=2, max=6, step=1, value=k_value),
        ui.button(name='predict_virus', label='Predict', primary=True)
    ])


    if 'predict_virus_tab' not in q.client.all_pages:
        q.client.all_pages.append('predict_virus_tab')

    #if q.client.promoter_data is not None:

    data_items = [ui.text_m(f'Loaded file "{q.client.virus_file_name}" has '
                            f'**{q.client.virus_data.shape[0]}** rows and **{q.client.virus_data.shape[1]}** features.\n\n'),
                  make_ui_table(q.client.virus_data, data_display_max_nrows)]
    q.page['virus_data_view'] = ui.form_card(box='4 2 9 4', items=data_items)

    if 'virus_data_view' not in q.client.all_pages:
        q.client.all_pages.append('virus_data_view')

    if q.client.virus_predictions is not None:

        download_data_text = '''=
Inference complete! Click [here]({{predictions}}) to download the predictions! Look below for visualization of virus composition and top kmers!
'''
        q.page['download_virus_predictions'] = ui.markdown_card(
                box='4 6 9 1',
                title='',
                content=download_data_text,
                data=dict(predictions=q.client.virus_predictions_path)
            )
        q.client.virus_data_predictions_seq_length= plot_virus_percent(q.client.virus_predictions)

        q.page['plot_virus'] = ui.frame_card(
            box='4 7 4 6',
            title='How many virus/non-virus are in the dataset',
            content=q.client.virus_data_predictions_seq_length
        )

        if 'download_virus_predictions' not in q.client.all_pages:
            q.client.all_pages.append('download_virus_predictions')

        if 'plot_virus' not in q.client.all_pages:
            q.client.all_pages.append('plot_virus')
        print(f'virus topk {q.client.virus_topk}')
        q.client.promoter_topk_plot=plot_top_kmers(q.client.top_virus_kmers, q.client.top_virus_kmer_counts, q.client.virus_topk)



        q.page['plot_virus_kmers'] = ui.frame_card(
            box='8 7 5 6',
            title='Top virus kmers extracted from attention weights',
            content=q.client.promoter_topk_plot
        )


        if 'plot_virus_kmers' not in q.client.all_pages:
            q.client.all_pages.append('plot_virus_kmers')


async def virus_model_predict(q):

    # Loading models
    start_time_main = time.time()

    if q.client.virus_models_loaded is None:
        q.page['progress'] = ui.form_card(box='4 6 9 1',
                                          items=[ui.progress(label="Loading models...")])
        await q.page.save()
        del q.page['progress']
        q.client.virus_inference=Virus_Inference.Virus_Inference()
        q.client.virus_inference.load_models('Virus_Inference/best_weights')
        q.client.virus_models_loaded=True

    # Getting predictions
    #if q.client.predictions is None:
    minimum_required_features = ["sequence"]
    #await progress_page(q,message="Getting predictions...")

    q.page['progress'] = ui.form_card(box='4 6 9 1',
                                      items=[ui.progress(label="Getting predictions...")])
    await q.page.save()
    del q.page['progress']

    q.client.virus_predictions, q.client.top_virus_kmers, q.client.top_virus_kmer_counts = \
    q.client.virus_inference.predict(q.client.virus_data.loc[:, minimum_required_features])
    q.client.virus_predictions.to_csv('temp/virus_predictions.csv', index=False)
    q.client.virus_predictions_path, = await q.site.upload(['temp/virus_predictions.csv'])

    elapsed_time = time.time() - start_time_main
    print(f"minutes passed for getting predictions: {round(elapsed_time/60, 2)}")




async def home(q):
    home_text = f'''

The Nucleic Transformer models are deep learning models developed to study and understand DNA/RNA usings public available datasets. You can check out [the paper on bioarxiv](https://www.biorxiv.org/content/10.1101/2021.01.28.428629v1)
and [open-sourced code on github](https://github.com/Shujun-He/Nucleic-Transformer). The model archiecture is simple but effective, outperforming previous results in DNA promoters/virus classification; additionally,
we used it to to place 7th in the [OpenVaccine challenge](https://www.kaggle.com/c/stanford-covid-vaccine).

Throughout the app, you will be able to use pretrained Nucleic Transformer models to classify DNA promoters/virus and quantify RNA degradation. Also,
Nucleic Transformer is very interpretible so you will be able to visualize the attention of the neural networks and understand what the neural network is looking at/for while making predictions.
This also allows extraction of k-mer promoter/viral motifs and informed decision making when using these neural networks.

![Plot]({q.app.home_image_1_url})

 '''
    q.page['home'] = ui.form_card(box='1 2 12 10',
        items=[
               ui.text_m(home_text)
               ])

async def display_data_parameters_page(q,random_sample_disabled=True):

    # Add navigation buttons
    print(q.args["#"])

    if q.client.activetab == "home":
            await delete_pages(q,keep_nav=True)
            await home(q)
    elif q.client.activetab == "promoter_prediction":
            await delete_pages(q,keep_nav=True)
            await predict_promoter_tab(q)
            print('line 546')

            #await predict_virus_tab(q)
    elif q.client.activetab == "enhancer_prediction":
            await delete_pages(q,keep_nav=True)
            await predict_enhancer_tab(q)
            print('line 546')


    elif q.client.activetab == "virusprediction":
            await delete_pages(q,keep_nav=True)
            await predict_virus_tab(q)
    elif q.client.activetab == "rnaprediction":
            await delete_pages(q,keep_nav=True)
            await predict_rna_tab(q)






async def main(q: Q):
    # Upload the logo & images
    if q.client.all_pages is None:
        q.client.all_pages = []

    if not q.app.logo_url:
        q.app.logo_url, = await q.site.upload([logo_file])

    if not q.app.home_image_1:
        q.app.home_image_1_url, = await q.site.upload([image_1_path])

    if not q.app.home_image_2:
        q.app.home_image_2_url, = await q.site.upload([image_2_path])

    if not q.app.home_image_3:
        q.app.home_image_3_url, = await q.site.upload([image_3_path])

    if not q.app.home_image_4:
        q.app.home_image_4_url, = await q.site.upload([image_4_path])


    if q.args.promoter_user_files:
        #await delete_pages(q)
        try:
            print('location: upload data')
            # Make the file available locally and store file path in client context
            q.client.local_path = await q.site.download(q.args.promoter_user_files[0], '.')
            q.client.link_to_file = q.args.promoter_user_files[0]

            if q.client.link_to_file.endswith('.json'):
                q.client.promoter_data = pd.read_json(q.client.local_path, lines=True)
            elif q.client.link_to_file.endswith('.csv'):
                q.client.promoter_data = pd.read_csv(q.client.local_path)
            q.client.fs_columns = list(q.client.promoter_data.columns.values.tolist())
            q.client.file_name = os.path.split(q.client.local_path)[1]
            q.client.target_columns = [col for col in target_columns if col in q.client.train.columns]
            q.client.promoter_data["sequence_length"] = q.client.promoter_data["sequence"].apply(lambda seq: len(seq))
            print(f"data shape: {q.client.promoter_data.shape}")

            data_items = [ui.text_m(f'Loaded file "{q.client.file_name}" has '
                                    f'**{q.client.promoter_data.shape[0]}** rows and **{q.client.train.shape[1]}** features.\n\n'),
                          make_ui_table(q.client.promoter_data, data_display_max_nrows)]
            q.page['promoter_data_view'] = ui.form_card(box='4 2 9 4', items=data_items)

            if 'promoter_data_view' not in q.client.all_pages:
                q.client.all_pages.append('promoter_data_view')

            if q.client.promoter_predictions is not None:
                q.client.promoter_predictions=None
                del q.page['plot_promoters']
                del q.page['plot_promoter_kmers']
                del q.page['download_promoter_predictions']

        #except Exception as e: print(e)
        except:
            await delete_pages(q, keep_nav=True)
            q.page['error'] = ui.form_card(box='1 2 6 2',
                items=display_error(error_type="upload"))

    elif q.args.enhancer_user_files:
        #await delete_pages(q)
        try:
            print('location: upload data')
            # Make the file available locally and store file path in client context
            q.client.local_path = await q.site.download(q.args.enhancer_user_files[0], '.')
            q.client.link_to_file = q.args.enhancer_user_files[0]

            q.client.enhancer_file=q.args.enhancer_user_files[0].split('/')[-1]

            if q.client.link_to_file.endswith('.json'):
                q.client.enhancer_data = pd.read_json(q.client.local_path, lines=True)
            elif q.client.link_to_file.endswith('.csv'):
                q.client.enhancer_data = pd.read_csv(q.client.local_path)
            q.client.fs_columns = list(q.client.enhancer_data.columns.values.tolist())
            q.client.enhancer_data["sequence_length"] = q.client.enhancer_data["sequence"].apply(lambda seq: len(seq))
            print(f"data shape: {q.client.enhancer_data.shape}")

            data_items = [ui.text_m(f'Loaded file "{q.client.enhancer_file}" has '
                                    f'**{q.client.enhancer_data.shape[0]}** rows and **{q.client.train.shape[1]}** features.\n\n'),
                          make_ui_table(q.client.enhancer_data, data_display_max_nrows)]
            q.page['enhancer_data_view'] = ui.form_card(box='4 2 9 4', items=data_items)
            #
            # if 'enhancer_data_view' not in q.client.all_pages:
            #     q.client.all_pages.append('enhancer_data_view')
            #
            # if q.client.enhancer_predictions is not None:
            #     q.client.enhancer_predictions=None
            #     del q.page['plot_enhancers']
            #     del q.page['plot_enhancer_kmers']
            #     del q.page['download_enhancer_predictions']

        #except Exception as e: print(e)
        except:
            await delete_pages(q, keep_nav=True)
            q.page['error'] = ui.form_card(box='1 2 6 2',
                items=display_error(error_type="upload"))


    elif q.client.local_path is None:
        #q.client.local_path = "stanford-covid-vaccine/train_small.json"
        q.client.local_path = "Promoter_Inference/promoter_sample.csv"
        #q.client.train = pd.read_json(q.client.local_path, lines=True)
        q.client.train = pd.read_csv(q.client.local_path)
        q.client.fs_columns = list(q.client.train.columns.values.tolist())
        q.client.file_name = os.path.split(q.client.local_path)[1]
        q.client.activetab = "home"
        q.client.target_columns = [col for col in target_columns if col in q.client.train.columns]
        q.client.train["sequence_length"] = q.client.train["sequence"].apply(lambda seq: len(seq))


        await display_data_parameters_page(q)
        await display_nav(q)
        #await display_file_upload(q)

    elif q.args.virus_user_files:
        #await delete_pages(q)
        try:
            print('location: upload data')
            # Make the file available locally and store file path in client context
            q.client.local_virus_file_path = await q.site.download(q.args.virus_user_files[0], '.')


            if q.client.local_virus_file_path.endswith('.json'):
                q.client.virus_data = pd.read_json(q.client.local_virus_file_path, lines=True)
            elif q.client.local_virus_file_path.endswith('.csv'):
                q.client.virus_data = pd.read_csv(q.client.local_virus_file_path)
            q.client.fs_columns = list(q.client.virus_data.columns.values.tolist())
            q.client.virus_file_name = os.path.split(q.client.local_virus_file_path)[1]
            q.client.target_columns = [col for col in target_columns if col in q.client.train.columns]
            q.client.virus_data["sequence_length"] = q.client.virus_data["sequence"].apply(lambda seq: len(seq))
            print(f"data shape: {q.client.virus_data.shape}")

            data_items = [ui.text_m(f'Loaded file "{q.client.virus_file_name}" has '
                                    f'**{q.client.virus_data.shape[0]}** rows and **{q.client.train.shape[1]}** features.\n\n'),
                          make_ui_table(q.client.virus_data, data_display_max_nrows)]
            q.page['virus_data_view'] = ui.form_card(box='4 2 9 4', items=data_items)

        # if 'virus_data_view' not in q.client.all_pages:
        #     q.client.all_pages.append('virus_data_view')

        # if q.client.virus_predictions is not None:
        #     q.client.virus_predictions=None
        #     del q.page['plot_virus']
        #     del q.page['plot_virus_kmers']
        #     del q.page['download_virus_predictions']

        #except Exception as e: print(e)
        except:
            await delete_pages(q, keep_nav=True)
            q.page['error'] = ui.form_card(box='1 2 6 2',
                items=display_error(error_type="upload"))

    elif q.client.local_virus_file_path is None:
        #q.client.local_path = "stanford-covid-vaccine/train_small.json"
        q.client.local_virus_file_path = "Virus_Inference/virus_sample.csv"
        #q.client.train = pd.read_json(q.client.local_path, lines=True)
        q.client.virus_data = pd.read_csv(q.client.local_virus_file_path)
        q.client.fs_columns = list(q.client.train.columns.values.tolist())
        q.client.virus_file_name = os.path.split(q.client.local_virus_file_path)[1]
        q.client.activetab = "home"
        q.client.target_columns = [col for col in target_columns if col in q.client.train.columns]
        q.client.virus_data["sequence_length"] = q.client.virus_data["sequence"].apply(lambda seq: len(seq))


        await display_data_parameters_page(q)
        await display_nav(q)
        #await display_file_upload(q)








        #This part is for controlling active tab page.
    elif q.args.predict_promoter:
        q.client.promoter_topk = q.args.promoter_topk
        await promoter_model_predict(q)
        print('line 815')
        await predict_promoter_tab(q)

    elif q.args.predict_enhancer:
        q.client.enhancer_topk = q.args.enhancer_topk
        await enhancer_model_predict(q)
        print('line 815')
        await predict_enhancer_tab(q)


    elif q.args.predict_rna:
        #q.client.virus_topk = q.args.virus_topk
        await rna_model_predict(q)
        print('line 815')
        await predict_rna_tab(q)

    elif q.args.predict_virus:
        q.client.virus_topk = q.args.virus_topk
        await virus_model_predict(q)
        print('line 815')
        await predict_virus_tab(q)

    elif q.args["#"]:
        q.client.activetab = q.args["#"]
        await display_data_parameters_page(q)

    else:
        print("location: else")
        await display_data_parameters_page(q)

    await q.page.save()


if __name__ == '__main__':
    listen('/nucleictransformer', main)
