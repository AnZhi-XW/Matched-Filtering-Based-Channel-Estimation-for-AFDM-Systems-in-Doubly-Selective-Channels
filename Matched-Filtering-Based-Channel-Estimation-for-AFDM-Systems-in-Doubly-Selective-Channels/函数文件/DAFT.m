function 	X	=	DAFT(x,M, c1, c2 )
%		x:              时域输入信号.
%		M:              子载波数.
%		c1 and c2:      DAFT参数.
%		X:              DAFT域输出信号.

num_Row		=	size( x, 1 );
num_Col		=	size( x, 2 );

if	nargin	<=	1
    M		=	num_Row;
end

if	M >= 1

    if	M == num_Row
        temp_Signal		=	x;
    elseif	M > num_Row
        temp_Signal		=	zeros( M, num_Col, 'like', x );
        temp_Signal( 1 : num_Row, : )	=	x;
    else
        temp_Signal	=	x( 1 : M, : );
        str_Warning_Msg	=	strcat( mfilename, ': The number of DAFT points is smaller than the number of rows of the input_signal' );
        warning( str_Warning_Msg );
    end


    chirp_Index	=	( 0 : M - 1 ).';
    Lamda1 = exp( -1i*2*pi*c1 * chirp_Index.^2 );
    Lamda2 = exp( -1i*2*pi*c2 * chirp_Index.^2 );

    x_temp1=repmat(Lamda1, 1, num_Col).*x;
    x_temp2=fft(x_temp1)./sqrt(M);
    x_temp3=repmat(Lamda2, 1, num_Col).*x_temp2;
    X=x_temp3;

end

end
