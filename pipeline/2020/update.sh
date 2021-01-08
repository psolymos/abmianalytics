# copy results
sudo cp -R /media/samba/abmisc/reports_staging/2020 /media/samba/abmisc/reports

# delete dist folder from staging dir
sudo rm -rf /media/samba/abmisc/reports_staging/dist

# delete website from root html
cd /var/www/html
sudo rm -rf _nuxt lichens mosses vplants mites birds mammals habitats

# copy over dist contents
sudo cp -R /media/samba/abmisc/reports_staging/dist/* /var/www/html
