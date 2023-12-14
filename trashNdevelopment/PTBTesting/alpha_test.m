grey = 128;
ms=100;
transLayer=2;
[x,y]=meshgrid(-ms:ms, -ms:ms);
maskblob=uint8(ones(2*ms+1, 2*ms+1, transLayer) * grey);
  
xsd=ms/2.0;
ysd=ms/2.0;
maskblob(:,:,transLayer)=uint8(round(255 - exp(-((x/xsd).^2)-((y/ysd).^2))*255));