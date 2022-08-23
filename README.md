# DAISGram
Repository for Daisgram project A.A. 2020-2021
## Library - BMP
To test the proper functioning of the library (libbmp) you need to
run the command `make testbmp` and then`. / Test_bmp`.

Image | Brighten (+20) | Brighten (+100) | Grayscale
------------ | ------------- | ------------- |-------------
![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/dais.bmp) | ![DAIS+20](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_brighten_20.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_brighten_100.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_gray.bmp)

Image | Smooth (h=3) | Smooth (h=5) | Smooth (h=7)
------------ | ------------- | ------------- | -------------
![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/dais.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_smooth_3.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_smooth_5.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_smooth_7.bmp)

Image | Sharp | Edge | Warhol
------------ | ------------- | ------------- | -------------
![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/dais.bmp) | ![DAIS+20](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_sharp.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_edge.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/dais_warhol.bmp) |

Image A | Image B | Blend alpha=0 | alpha=0.25 | alpha=0.5 | alpha=0.75 | alpha=1.00
------------ | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/blend/blend_a.bmp) | ![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/blend/blend_b.bmp) | ![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/results/blend/blend_0.00.bmp) | ![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/results/blend/blend_0.25.bmp) | ![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/results/blend/blend_0.50.bmp) | ![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/results/blend/blend_0.75.bmp) | ![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/results/blend/blend_1.00.bmp)

## Green Screen
Close-Up Piano | Background | GreenScreen
------------ | ------------- | -------------
![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/greenscreen/gs_2.bmp) | ![DAIS+20](https://github.com/RiccardoSale/DAISgram/blob/main/images/greenscreen/gs_2_bkg.bmp) | ![DAIS+100](https://github.com/RiccardoSale/DAISgram/blob/main/results/greenscreen/dais_matrix.bmp)

 ## Histogram equalization
 Original Image | Equalization
------------ | -------------
![DAIS](https://github.com/RiccardoSale/DAISgram/blob/main/images/equalize/fullmoon.bmp) | ![DAIS+20](https://github.com/RiccardoSale/DAISgram/blob/main/results/fullmoon_equalize.bmp)
