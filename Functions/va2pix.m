function va_in_pixels = va2pix(distance,px_per_cm)  
  # converts a unit visual angle to number of pixels on the screen
  # 2023 Kevin Thai, Lauri Nurminen
  va = tan(pi/180) * distance;
  va_in_pixels = va * px_per_cm;
end