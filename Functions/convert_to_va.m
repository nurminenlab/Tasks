function va_in_pixels = convert_to_va(distance)
  cm_to_pixels = 37.8;
  va = tan(pi/180) * 2 * distance;
  va_in_pixels = va * cm_to_pixels;
end