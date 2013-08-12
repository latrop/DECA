echo ""
echo "******** IRAF INSTALLATION SCRIPT FOR LINUX ********"

os=$(uname -s)
case $os in
  Linux) 
    distro=$(lsb_release -si)
    case $distro in
      Ubuntu) 
        echo "Ubuntu, great!"
        ;;
      *) 
        echo "Sorry, Linux distribution '$distro' is not supported"
        exit 1
        ;;
    esac
    ;;
esac


OS=$(uname -m)
case $OS in
  x86_64) 
    echo "Downloading IRAF for x86_64..."
    #wget ftp://iraf.noao.edu/iraf/v216/PCIX/iraf.lnux.x86_64.tar.gz	# for 64
    mv iraf.lnux.x86_64.tar.gz iraf.lnux.tar.gz
    ;;
  *) 
    echo "Downloading IRAF for 32 architecture..."
    #wget ftp://iraf.noao.edu/iraf/v216/PCIX/iraf.lnux.x86.tar.gz	# for 32
    mv iraf.lnux.x86.tar.gz iraf.lnux.tar.gz
    ;;
esac

echo "Downloading x11iraf-v2.0BETA..."
#wget http://iraf.noao.edu/x11iraf/x11iraf-v2.0BETA-bin.linux.tar.gz


sudo -s
mkdir /iraf
mkdir /iraf/iraf
mkdir /iraf/x11iraf

mv x11iraf-v2.0BEAT-bin.linux.tar.gz /iraf/x11iraf/
mv iraf.lnux.tar.gz /iraf/iraf/

cd /iraf/iraf
tar zxvf iraf.lnux.tar.gz

cd /iraf/x11iraf/
tar zxvf x11iraf-v2.0BETA-bin.linux.tar.gz

apt-get install tcsh

# Here you need to enter a new password for IRAF user and confirm it. Leave all other questions without answers (ENTER)! When they ask if the information is correct, type Y.
adduser iraf --home /iraf/iraf/local --shell /bin/csh --ingroup sudo

chown iraf -R /iraf
su iraf
sudo -s
setenv iraf /iraf/iraf
cd /iraf/iraf/unix/hlib
source irafuser.csh

case $OS in
  x86_64)
    echo "Installing 32 libraries..."
    apt-get install ia32-libs
    ;;
  *) 
    echo "There is no need to install 32 libraries!"
    ;;
esac

echo "X11IRAF INSTALLATION!"
cd /iraf/x11iraf/
./install
# Answer the following questions:
# Type /usr/local/bin
# Yes
# Enter
# Enter
# Enter
# ENTER
# Enter


echo "IRAF INSTALLATION!"
cd /iraf/iraf/unix/hlib
#./install -n
# Enter
# Enter
# Enter
# Enter
# Enter
# Enter
# Configure IRAF Networking on this machine: no !!!!!!!!!!!!!!!!!!!!!!!
# Enter
# Enter
# Enter
# Enter
./install
# Enter
# Enter
# Enter
# Enter
# Enter
# Enter
# Configure IRAF Networking on this machine: no !!!!!!!!!!!!!!!!!!!!!!!
# Enter
# Enter
# Enter
# Enter

echo "STSDAS INSTALLATION!"
cd ~
mkdir IRAF
mkdir IRAF/Packages
mkdir IRAF/Packages/STSDAS
cd IRAF/Packages/STSDAS/
wget http://stsdas.stsci.edu/download/release_2013-03/stsci_iraf-3.16.redhat.tar.gz
tar zxvf stsci_iraf-3.16.redhat.tar.gz
cd stsci_iraf-3.16/
su iraf
sudo -s
./install_helper
echo "Follow the instructions you've just received! Copy the upper lines and paste them to the indicated IRAF file! "






echo "******** THE INSTALLATION IS DONE! ********"
