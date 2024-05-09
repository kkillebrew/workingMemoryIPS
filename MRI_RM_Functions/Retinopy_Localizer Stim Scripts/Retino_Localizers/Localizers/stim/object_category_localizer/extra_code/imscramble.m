%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Scramble                      %
%%              Chris Longmore Sep 04               %
%%               Univeristy of York                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reads in an image and scrambles it up in a series
%of rectangles
%Usage:
%  scramble('in_name','out_name', x, y) where:
%  in_name = name of input file
%  out_name = name of output file
%  x = number of blocks to create horizontally
%  y = number of blocks to create vertially
% NOTES:
% x must be a equal divisor of the width of the input image
% y must be a equal divisor of the height of the input image
% Filenames must include an extension (e.g. 'fname.bmp' will work but
% 'fname' will not

function imscramble(fname, oname, y, x)

[pathstr,name,ext,versn] = fileparts(fname);                        %Get the properties of the input filename

if strcmp(ext,'')                                                   %Check that the extension is present and if not...
    error ('Please specify the input filename WITH the extension')  %Inform the user and quit
end

[pathstr,name,ext,versn] = fileparts(oname);                        %Get the properties of the output filename

if strcmp(ext,'')                                                   %Check that the extension is present and if not...
    error ('Please specify the output filename WITH the extension') %Inform the user and quit
end

imginfo = imfinfo(fname);

yres = getfield(imginfo,'Width');                                   %Get width of image
xres = getfield(imginfo,'Height');                                  %Get height of image

if getfield(imginfo,'BitDepth') == 24                               %Get depth of image
    truecolour = true;
    scramimg = uint8(zeros(xres,yres,3));                           %Create the empty matrix for the scrambled image
else
    truecolour = false;
    scramimg = uint8(zeros(xres,yres));                             %Create the empty matrix for the scrambled image
end

if truecolour                                                       %If the image is 24 bit then read image into x-y-3 matrix
    img = imread(fname);
else
    [img,map] = imread(fname);                                      %Else read into x-y matrix with a colourmap
end

fprintf ('Input image is %g pixels by %g pixels, %g-bit\n',yres, xres, getfield(imginfo,'BitDepth'))      %Tell the user the imput image dimenstions

if rem(xres,x) ~= 0                                                 %Check the height is divisable by the number of tiles requested
    error ('Height of image is not divisable by the number of tiles requested')
end

if rem(yres,y) ~= 0                                                 %Check the width is divisable by the number of tile requested
    error ('Width of image is not divisable by the number of tiles requested')
end

distx = xres / x;                                                   %Get how many pixels each tile is high
disty = yres / y;                                                   %Get how many pixels each tile is wide
    
z = x * y;                                                          %Figure out how many tiles there will be altogether

%The logic behind this program is to divide the original image up into a
%series of tiles.  Each tile is numbered 1 to z (which equals the total)
%number of tiles needed to divide the image up.  We first generate an
%array containing the numbers of the tile and then shuffle them into a
%random order.  This new order is then used to construct the new image.
%For example, if we request an image consisting of 2 x 2 tiles then we have
%4 tiles; 1, 2, 3 and 4.  This order is shuffled (to say; 4, 2, 1 and 3).
%In the new images, tile 1 is tile 4 of the old image, tile 2 is tile 2
%from the old image, tile 3 is tile 1 and tile 4 is tile 3.  Finally, the
%image is written to a new file.

for i=1:z                                                           %Create references for z tiles
    order(i) = i;
end

order = Shuffle(order);                                             %Shuffle the order of the tiles
ok = false;                                                         %Needs for loop below
a = 1;                                                              %Starting Y value for the new image
b = 1;                                                              %Starting X value for the new image
row = 0;                                                            %Which row we are currently writing (always real row - 1)

for loop = 1:z                                                      %Start main loop
    count = 1;                                                      %Counter to determine which column our target tile in the original image is in
    tmp = order(loop);                                              
    while ~ok                                                       %Determine which column the target tile is on
        if tmp <= x
            ok = true;
        else
            count = count + 1;
            tmp = tmp - x;
        end
    end
    xpos = order(loop) - (x * (count - 1));                         %Calculate target column
    ypos = count;                                                   %Automatically get the target row
    ok = false;
    for i=(distx*(xpos-1))+1:(distx*xpos)                           %Move to the point of our target tile in the original image
        for j=(disty*(ypos-1))+1:(disty*ypos)
            if truecolour
                scramimg(a,b,1) = img(i,j,1);                           %And write it into the location in the new image
                scramimg(a,b,2) = img(i,j,2);
                scramimg(a,b,3) = img(i,j,3);
            else
                scramimg(a,b) = img(i,j);                           %For colourmapped images
            end
            b = b + 1;
        end
        b = (disty * row) + 1;
        a = a + 1;
    end
    if a > xres                                                     %Reach the end of the row in the new image so move to the next one
        a = 1;
        row = row + 1;
        b = (disty * row)+1;
    end
end

fprintf ('Writing output file\n')
if truecolour
    imwrite(scramimg,oname);                                        %Write out scrambled image
else
    imwrite(scramimg,map,oname)
end

fprintf ('Scramble successful\n')