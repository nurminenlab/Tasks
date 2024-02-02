## Copyright (C) 2023 nurminenlab
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## Author: Lauri Nurminen <lonurminen@uh.edu>
## Created: 2023-11-29
function Trans_mat = manual_translate(trial_records,fixation_rect,window_pointer,texture_pointer,Trans_mat)

  [center_x,center_y] = WindowCenter(window_pointer);
  small_rect  = [0 0 5 5];
  small_rects = ones(4,size(trial_records(end).gaze_position,2));
  for i = 1:size(trial_records(end).gaze_position,2)
      small_rects(:,i) = CenterRectOnPoint(small_rect,trial_records(end).gaze_position(1,i),trial_records(end).gaze_position(2,i))';
  endfor 

  Screen('DrawTexture', window_pointer, texture_pointer, [], fixation_rect, 0);
  Screen('FillOval', window_pointer, [255, 0, 0], small_rects);
  Screen('Flip', window_pointer);
  [clicks,x,y,whichButton] = GetClicks(window_pointer, 0.1);    
  x = x - center_x;
  y = y - center_y;
  
  Trans_mat = [-x,-y]';  
  for i = 1:size(trial_records(end).gaze_position,2)
      [x, y] = RectCenter(small_rects(:,i));
      small_rects(:,i) = CenterRectOnPoint(small_rects(:,i),x+Trans_mat(1),y+Trans_mat(2))';      
  endfor
  
  Screen('DrawTexture', window_pointer, texture_pointer, [], fixation_rect, 0);
  Screen('FillOval', window_pointer, [0, 0, 255], small_rects);
  Screen('Flip', window_pointer);
  KbWait();
endfunction
