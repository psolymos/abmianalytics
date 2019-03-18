# SC drive setup

## Login

```
ssh <user>@<IP_address>
```

```
sudo apt-get update
#sudo apt-get dist-upgrade
sudo apt-get upgrade

sudo -i
echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt-get update

sudo apt-get install r-base
sudo apt-get install r-base-dev

sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get install -y \
    pandoc \
    pandoc-citeproc \
    build-essential \
    libssl-dev \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxml2-dev \
    libxt-dev \
    libv8-dev \
    liblwgeom-dev \
    git \
    jags \
    ruby \
    ruby-dev \
    make \
    gcc \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    golang-go \
    nodejs \
    npm

## build jekyll
sudo gem install jekyll bundler

## set up webhooks
export GOPATH=$HOME/go
go get github.com/adnanh/webhook

## Apache (OpenCPU gets this done)
sudo apt install apache2
sudo ufw app list
sudo ufw allow 'Apache'
sudo ufw status
sudo systemctl status apache2

## OCPU
sudo add-apt-repository -y ppa:opencpu/opencpu-2.1
sudo apt-get update 
sudo apt-get upgrade
# Installs OpenCPU server
sudo apt-get install -y opencpu-server


## Rstudio server
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb
sudo gdebi rstudio-server-1.1.463-amd64.deb
sudo ufw allow 8787

## shiny-server
sudo su - \
-c "R -e \"install.packages('shiny', repos='https://cran.rstudio.com/')\""
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb
sudo gdebi shiny-server-1.5.9.923-amd64.deb
sudo ufw allow 3838



sudo apt-get install docker
sudo apt-get install docker-compose

```

## Samba

https://linuxize.com/post/how-to-install-and-configure-samba-on-ubuntu-18-04/

```
sudo apt update
sudo apt install samba
sudo systemctl status nmbd
sudo ufw allow 'Samba'

sudo mkdir /samba
sudo chgrp sambashare /samba
## add admin user group 'sadmin'
sudo useradd -M -d /samba/users -s /usr/sbin/nologin -G sambashare sadmin
sudo smbpasswd -a sadmin
sudo smbpasswd -e sadmin
sudo mkdir /samba/users
sudo chown sadmin:sambashare /samba/users
sudo chmod 2770 /samba/users
```

Add this by `sudo vi /etc/samba/smb.conf`:

```
[users]
    path = /samba/users
    browseable = yes
    read only = no
    force create mode = 0660
    force directory mode = 2770
    valid users = @sambashare @sadmin
```

Then restart: `sudo systemctl restart nmbd`