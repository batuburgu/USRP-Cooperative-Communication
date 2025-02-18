function [error]= error_calculation(complex_number)

%not sure how it works but the a and b values need to be changed based on your
%constellation phase, wanted to try writing cos(pi/4) for our case however, got
%chickened out too much to try, leave as it is just download it and do not
%touch this function again :D


    if(real(complex_number)>0)
        a=sqrt(2)/2;
    else
        a=-sqrt(2)/2;
    end


    if(imag(complex_number)>0)
        b=sqrt(2)/2;
    else
        b=-sqrt(2)/2;
    end

    error=a*imag(complex_number)-b*real(complex_number);
end
