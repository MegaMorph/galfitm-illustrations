cd models
python make_model_feedmes.py
python make_images.py
cd ../fits
python make_fit_feedmes.py
python run_fits.py
python make_aperture_feedmes.py
python run_fits.py '*[36]'
cd ..
python plot_illustrations.py
