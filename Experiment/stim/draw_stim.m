% shared stimulus drawing code for the sample and probe periods.
% this way, we don't have to replicate this bit of code in two places,
% especially given the switch statement to handle different stimtypes.

% N.B. be sure to set this.stim_phase appropriately before calling
%      draw_stim

% determine if target defined for sample or probe phase
this.targ_feature = this.(['targ_feature_' this.stim_phase]); % for easy reference below

switch trials.stimtype{it,1}
    case 'orientation' % oriented bars
        % start with template for a vertically oriented bar
        vert_rect = stim.orientation.rect;
        % draw targets
        for j = 1:size(this.targ_feature,1)
            cur_rect = RotatePointsPTB(vert_rect,this.targ_feature(j,1),[xc yc]); % rotate bar based on targ_feature
            r = CenterRectOnPoint(cur_rect,this.targ_pos_pix(j,1),this.targ_pos_pix(j,2)); % translate to appropriate coordinate
            Screen('FillRect',w,params.stim.orientation.targ_color,r);
        end
        % draw distractors
        for j = 1:size(this.dist_feature,1)
            cur_rect = RotatePointsPTB(vert_rect,this.dist_feature(j,1),[xc yc]); % rotate bar based on dist_feature
            r = CenterRectOnPoint(cur_rect,this.dist_pos_pix(j,1),this.dist_pos_pix(j,2)); % translate to appropriate coordinate
            Screen('FillRect',w,params.stim.orientation.dist_color,r);
        end
        
    case 'color' % colored circles/squares
        % start with template for any shape
        temp_rect = stim.color.rect;
        % draw targets
        for j = 1:size(this.targ_feature,1)
            switch params.stim.color.targ_shape
                case 'circle'
                    tmp_draw_cmd = 'FillOval';
                case 'square'
                    tmp_draw_cmd = 'FillRect';
                otherwise
                    error('don''t know which draw command to use for a targ_shape of %s',params.stim.color.targ_shape)
            end
            r = CenterRectOnPoint(temp_rect,this.targ_pos_pix(j,1),this.targ_pos_pix(j,2)); % translate to appropriate coordinate
            Screen(tmp_draw_cmd,w,this.targ_feature(j,:),r);
        end
        % draw distractors
        for j = 1:size(this.dist_feature,1)
            switch params.stim.color.dist_shape
                case 'circle'
                    tmp_draw_cmd = 'FillOval';
                case 'square'
                    tmp_draw_cmd = 'FillRect';
                otherwise
                    error('don''t know which draw command to use for a dist_shape of %s',params.stim.color.dist_shape)
            end
            r = CenterRectOnPoint(temp_rect,this.dist_pos_pix(j,1),this.dist_pos_pix(j,2)); % translate to appropriate coordinate
            Screen(tmp_draw_cmd,w,this.dist_feature(j,:),r);
        end
        
    case 'squares' % squares and diaomnds
        % start with point list for a square
        sq_points = stim.squares.points;
        % draw targets
        for j = 1:size(this.targ_feature,1)
            cur_rect = RotatePointsPTB(sq_points,this.targ_feature(j,1),[xc yc]); % rotate bar based on targ_feature
            r = bsxfun(@plus,cur_rect,this.targ_pos_pix(j,:)-[xc yc]); % translate point list to appropriate coordinate
            Screen('FillPoly',w,params.stim.squares.targ_color,r,1);
        end
        % draw distractors
        for j = 1:size(this.dist_feature,1)
            cur_rect = RotatePointsPTB(sq_points,this.dist_feature(j,1),[xc yc]); % rotate bar based on dist_feature
            r = bsxfun(@plus,cur_rect,this.dist_pos_pix(j,:)-[xc yc]); % translate point list to appropriate coordinate
            Screen('FillPoly',w,params.stim.squares.dist_color,r,1);
        end

    case 'letters' % colored letters
        % set text size and font for stimuli
        old_textSize = Screen('TextSize',w,params.stim.letters.letter_size);
        old_textFont = Screen('TextFont',w,params.stim.letters.letter_font);
        % draw targets
        for j = 1:size(this.targ_feature,1)
            tmp_letter = this.targ_feature(j,1); % get this letter
            tmp_bounds = Screen('TextBounds',w,tmp_letter);
            tmp_text_width  = RectWidth(tmp_bounds);
            tmp_text_height = RectHeight(tmp_bounds);
            Screen('DrawText',w,tmp_letter,this.targ_pos_pix(j,1)-tmp_text_width/2,this.targ_pos_pix(j,2)-tmp_text_height/2,params.stim.letters.targ_color);
        end
        % draw distractors
        for j = 1:size(this.dist_feature,1)
            tmp_letter = this.dist_feature(j,1); % get this letter
            tmp_bounds = Screen('TextBounds',w,tmp_letter);
            tmp_text_width  = RectWidth(tmp_bounds);
            tmp_text_height = RectHeight(tmp_bounds);
            Screen('DrawText',w,tmp_letter,this.dist_pos_pix(j,1)-tmp_text_width/2,this.dist_pos_pix(j,2)-tmp_text_height/2,params.stim.letters.dist_color);
        end
        % reset text size and font
        Screen('TextSize',w,old_textSize);
        Screen('TextFont',w,old_textFont);
        
    otherwise
        error('sorry, i don''t know how to draw stimuli of stimtype %s',trials.stimtype{it,1})
end