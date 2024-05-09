clear

 [id, names]  = GetKeyboardIndices
%[id names]  = GetGamepadIndices

% idx = strcmp('fORP Interface',names)
% keyboard_index = id(idx)
keyboard_index= str2num(input('Press enter the keyboard index for button presses','s'))
trigger_index= str2num(input('Press enter the keyboard index for scanner triggers','s'))

datafile = input('Please Hit Enter','s');

repetitions=1;

ListenChar(2);
HideCursor;
[w,rect]=Screen('OpenWindow',0,[128 128 128]);
xc=rect(3)/2;
yc=rect(4)/2;

one=KbName('1!');
two=KbName('2@');
three=KbName('3#');
four=KbName('4$');
escape = KbName('escape');

redfont = [255 0 0];
blackfont = [0 0 0];
greenfont = [0 255 0];

condition_list=[];

conditions=[1 2 3 4];

num_trials=repetitions*length(conditions);

for i=1:repetitions
    condition_list=[condition_list;fullyfact([4])];
end



trial_order=randperm(num_trials);
% 
        WaitSecs(1);
        Screen('TextSize',w,24);
        text3='When given the cue, you will press each button on the keypad from left to right, one at a time';
        width=RectWidth(Screen('TextBounds',w,text3));
        Screen('DrawText',w,text3,xc-width/2,yc,[0 0 0]);
        Screen('TextSize',w,24);
        text3='Please wait for specific instructions to press each key';
        width=RectWidth(Screen('TextBounds',w,text3));
        Screen('DrawText',w,text3,xc-width/2,yc+50,[0 0 0]);
        Screen('TextSize',w,24);
        Screen('TextSize',w,24);
        text3='When ready, press any key to continue';
        width=RectWidth(Screen('TextBounds',w,text3));
        Screen('DrawText',w,text3,xc-width/2,yc+100,[0 0 0]);
        Screen('TextSize',w,24);
        Screen('Flip',w);      
        KbWait(keyboard_index);
        
        
        WaitSecs(.5);
        Screen('TextSize',w,24);
        width=RectWidth(Screen('TextBounds',w,'Please press the leftmost key on the keypad'));
        Screen('DrawText',w,'Please press the leftmost key on the keypad',xc-width/2,yc,[0 0 0]);
        Screen('Flip',w);
        
        while KbCheck(keyboard_index); end 
        
        while 1
            
            [ keyIsDown, seconds, keyCode ] = KbCheck(keyboard_index);
            
            if keyIsDown
                button1 = KbName(keyCode); 
               
                WaitSecs(.1);
                Screen('TextSize',w,24);
                width=RectWidth(Screen('TextBounds',w,sprintf('%s%s', 'You pressed ', button1)));
                Screen('DrawText',w,sprintf('%s%s', 'You pressed ', button1),xc-width/2,yc,[0 0 0]);
                Screen('Flip',w);
                
                button1 = KbName(keyCode);                 
                
                while KbCheck(keyboard_index); end
                break
            end
            
        end

        
        
        
        
        WaitSecs(.5);
        Screen('TextSize',w,24);
        width=RectWidth(Screen('TextBounds',w,'Please press the next key going right on the keypad'));
        Screen('DrawText',w,'Please press the next key going right on the keypad',xc-width/2,yc,[0 0 0]);
        Screen('Flip',w);
        while KbCheck(keyboard_index); end 
        
        while 1
            
            [ keyIsDown, seconds, keyCode ] = KbCheck(keyboard_index);
            
     
            if keyIsDown
                button2 = KbName(keyCode);
               
                WaitSecs(.5);
                Screen('TextSize',w,24);
                width=RectWidth(Screen('TextBounds',w,sprintf('%s%s', 'You pressed ', button2)));
                Screen('DrawText',w,sprintf('%s%s', 'You pressed ', button2),xc-width/2,yc,[0 0 0]);
                Screen('Flip',w);
                
                
%            
                while KbCheck(keyboard_index); end
                break
            end
            
        end
        
        WaitSecs(.5);
        Screen('TextSize',w,24);
        width=RectWidth(Screen('TextBounds',w,'Please press the next key after that'));
        Screen('DrawText',w,'Please press the next key going right on the keypad',xc-width/2,yc,[0 0 0]);
        Screen('Flip',w);
        while KbCheck(keyboard_index); end 
        
        while 1
            
            [ keyIsDown, seconds, keyCode ] = KbCheck(keyboard_index);           
            
            if keyIsDown
                
                button3 = KbName(keyCode);
                
                WaitSecs(.5);
                Screen('TextSize',w,24);
                width=RectWidth(Screen('TextBounds',w,sprintf('%s%s', 'You pressed ', button3)));
                Screen('DrawText',w,sprintf('%s%s', 'You pressed ', button3),xc-width/2,yc,[0 0 0]);
                Screen('Flip',w);
                
                
%                 
                while KbCheck(keyboard_index); end
                break
            end
            
        end
        
        WaitSecs(.5);
        Screen('TextSize',w,24);
        width=RectWidth(Screen('TextBounds',w,'Please press the last rightmost key on the keypad'));
        Screen('DrawText',w,'Please press the last rightmost key on the keypad',xc-width/2,yc,[0 0 0]);
        Screen('Flip',w);
        while KbCheck(keyboard_index); end 
        
        while 1
            
            [ keyIsDown, seconds, keyCode ] = KbCheck(keyboard_index);
            
            if keyIsDown
                button4 = KbName(keyCode);
          
                WaitSecs(.5);
                Screen('TextSize',w,24);
                width=RectWidth(Screen('TextBounds',w,sprintf('%s%s', 'You pressed ', button4)));
                Screen('DrawText',w,sprintf('%s%s', 'You pressed ', button4),xc-width/2,yc,[0 0 0]);
                Screen('Flip',w);
                
                
                 
                while KbCheck(keyboard_index); end
                break
            end
            
        end
        
        
%         WaitSecs(1);
%         Screen('TextSize',w,24);
%         text5='Please press the appropriate keys as will be indicted by the instructions';
%         width=RectWidth(Screen('TextBounds',w,text3));
%         Screen('DrawText',w,text5,xc-width/2,yc,[0 0 0]);
%         Screen('TextSize',w,24);
%         text5='When ready, press any key to continue';
%         width=RectWidth(Screen('TextBounds',w,text3));
%         Screen('DrawText',w,text5,xc-width/2,yc+50,[0 0 0]);
%         Screen('TextSize',w,24);     
%         KbWait;
        

for i=1:num_trials
    
    condition=condition_list(trial_order(i));
    moveon=0;
    switch condition
        case 1
            text= sprintf('%s%s', 'Please press the leftmost key again: ',  button1);
        case 2
            text= sprintf('%s%s', 'Please press the second key from the left: ',  button2);
        case 3
            text= sprintf('%s%s', 'Please press the third key from the left: ',  button3);
        case 4
            text= sprintf('%s%s', 'Please press the rightmost key: ',  button4);
    end
    
        
    
    while moveon==0;
        WaitSecs(.5)
        Screen('TextSize',w,24);
        width=RectWidth(Screen('TextBounds',w,text));
        Screen('DrawText',w,text,xc-width/2,yc,[0 0 0]);
        Screen('Flip',w);
        [keyisdown, secs, keycode] = KbCheck(keyboard_index);
        while 1
            [keyisdown, secs, keycode] = KbCheck(keyboard_index);
            switch condition
                case 1
                    if keycode(KbName(button1));
                        moveon=1;
                        break
                    elseif keyisdown==1;
                        break
                    end
                    
                case 2
                    if keycode(KbName(button2));
                        moveon=1;
                        break
                    elseif keyisdown==1;
                        break
                    end
                    
                case 3
                    if keycode(KbName(button3));
                        moveon=1;
                        break
                    elseif keyisdown==1;
                        break
                    end
                    
                case 4
                    if keycode(KbName(button4));
                        moveon=1;
                        break
                    elseif keyisdown==1;
                        break
                    end
            end
            if keycode(escape)
                Screen('CloseAll')
            end
        end
        if moveon~=1
            text2='I am sorry. Please try again.';
            font = redfont;
        else
            text2='Correct!';
            font = greenfont;
        end
        Screen('TextSize',w,24);
        
        width=RectWidth(Screen('TextBounds',w,text2));
        Screen('DrawText',w,text2,xc-width/2,yc,font);
        Screen('Flip',w);
        WaitSecs(1);
        
        
    end
    
end

KbWait(trigger_index);


ShowCursor;
ListenChar(0);
Screen('CloseAll');