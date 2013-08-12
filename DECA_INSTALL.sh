echo ""
echo "************************************************************"
echo "******** DECA INSTALLATION SCRIPT FOR UBUNTU/DEBIAN ********"

echo "Please, choose the packages to install:"
echo ""
echo "	0. DECA"
echo "	1. SciPy with NumPy, IPython, MatPlotLib"
echo "	2. PyFITS (pip-installation)"
echo "	3. fpdf (source)"
echo "	4. ds9 (executable)"
echo "	5. SExtractor (deb32/64)"
echo "	6. GALFIT (deb32/64)"
echo "	7. IRAF with STSDAS"
echo "	8. readAtlasImages from SDSS (source)"
#echo "	9. All of the above"
echo "	10. Exit"
echo ""







while true; do
    read -p "Enter the number of action?" act
    case $act in
	"0" | "9" )
		# 0. DECA INSTALLATION
		echo "The PATH willbe added to your .bachrc, so you can call DECA just by typing $ deca.py OR you can do it by hand."
		while true; do
		    read -p "Do you want to continue?" yn
		    case $yn in
			[Yy]* ) 
				path=$(pwd)
				cd ~
				echo "export PATH="$path":$""PATH" >> .bashrc
				cd $path;;
			[Nn]* )
				echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;

	"1" | "9" )
		# 1.SciPy with NumPy, IPython, MatPlotLib
		while true; do
		    read -p "Do you wish to install SciPy with NumPy, IPython, MatPlotLib?" yn
		    case $yn in
			[Yy]* )
				sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose
				break;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	"2" | "9" )
		# 2. PyFITS
		sudo apt-get install python-pip
		while true; do
		    read -p "Do you wish to install PyFITS?" yn
		    case $yn in
			[Yy]* ) pip install pyfits;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		done;;


	"3" | "9" )
		# 3. fpdf
		while true; do
		    read -p "Do you wish to install fpdf?" yn
		    case $yn in
			[Yy]* ) wget https://code.google.com/p/pyfpdf/downloads/detail?name=fpdf-1.7.zip
				unzip fpdf-1.7.zip
				cd fpdf-1.7
				sudo python setup.py install;;

			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	"4" | "9" )
		#4. ds9
		OS=$(uname -m)
		while true; do
		    read -p "Do you wish to install ds9?" yn
		    case $yn in
			[Yy]* ) 
				case $OS in
				  x86_64) 
				    echo "Downloading ds9 for x86_64..."
				    wget http://hea-www.harvard.edu/RD/ds9/download/linux64/ds9.linux64.7.2.tar.gz	# for 64
				    tar zxvf ds9.linux64.7.2.tar.gz
				    sudo mv ds9 /usr/local/bin
				    sudo chmod +x /usr/local/bin/ds9
				    sudo apt-get install libxss1
				    ;;
				  *) 
				    echo "Downloading ds9 for 32 architecture..."
				    wget http://hea-www.harvard.edu/RD/ds9/download/linux/ds9.linux.7.2.tar.gz	# for 32
				    tar zxvf ds9.linux.7.2.tar.gz
				    sudo mv ds9 /usr/local/bin
				    sudo chmod +x /usr/local/bin/ds9
				    sudo apt-get install libxss1
				    ;;
				esac;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	"5" | "9" )
		# 5. SExtractor
		while true; do
		    read -p "Do you wish to install SExtractor?" yn
		    case $yn in
			[Yy]* ) 
				case $OS in
				  x86_64) 
				    echo "Installing SExtractor..."
				    wget http://pkgs.org/ubuntu-12.04/ubuntu-universe-amd64/sextractor_2.4.4-1_amd64.deb/download/	# for 64
				    sudo apt-get update
				    sudo apt-get install sextractor
				    ;;
				  *) 
				    echo "Installing SExtractor..."
				    wget http://pkgs.org/ubuntu-12.04/ubuntu-universe-i386/sextractor_2.4.4-1_i386.deb/download/	# for 32
				    sudo apt-get update
				    sudo apt-get install sextractor
				    ;;
				esac;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	"6" | "9" )
		# 6. GALFIT
		while true; do
		    read -p "Do you wish to install GALFIT?" yn
		    case $yn in
			[Yy]* ) 
				case $OS in
				  x86_64) 
				    echo "Installing GALFIT..."
				    wget http://users.obs.carnegiescience.edu/peng/work/galfit/galfit3-debian64.tar.gz	# for 64
				    tar zxvf galfit3-debian64.tar.gz
				    sudo mv galfit /usr/local/bin
				    sudo chmod +x /usr/local/bin/galfit
				    ;;
				  *) 
				    echo "Installing GALFIT..."
				    wget http://users.obs.carnegiescience.edu/peng/work/galfit/galfit3-debian32.tar.gz	# for 32
				    tar zxvf galfit3-debian32.tar.gz
				    sudo mv galfit /usr/local/bin
				    sudo chmod +x /usr/local/bin/galfit
				    ;;
				esac;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	"7" | "9" )
		# 7. IRAF with STSDAS
		while true; do
		    read -p "Do you wish to install IRAF with STSDAS?" yn
		    case $yn in
			[Yy]* ) sh IRAF_INSTALL.sh;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	"8" | "9" )
		# 8. readAtlasImages
		while true; do
		    read -p "Do you wish to install readAtlasImages from SDSS?" yn
		    case $yn in
			[Yy]* ) 
				wget http://www.sdss3.org/binaries/readAtlasImages-v5_4_11.tar.gz
				tar -xvf readAtlasImages-v5_4_11.tar.gz
				cd readAtlasImages-v5_4_11
				make clean
				make
				sudo mv read_PSF /usr/local/bin
				;;
			[Nn]* ) echo "The installation is missed."
				exit 0;;
			* ) echo "Please answer yes or no.";;
		    esac
		    break
		    echo ""
		done;;


	10) exit 0;;
    esac
done









