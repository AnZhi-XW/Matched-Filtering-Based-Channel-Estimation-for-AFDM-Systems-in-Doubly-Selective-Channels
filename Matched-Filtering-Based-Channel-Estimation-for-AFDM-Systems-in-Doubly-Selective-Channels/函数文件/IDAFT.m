function 	x	=	IDAFT(X,M, c1, c2 )
%		X:              DAFT域输入信号
%		M:              子载波数
%		c1 and c2:      IDAFT参数
%		x:              时域输出信号

num_Row		=	size( X, 1 );
num_Col		=	size( X, 2 );

if	nargin	<=	1
    M		=	num_Row;
end

if	M >= 1

    if	M == num_Row
        temp_Signal		=	X;
    elseif	M > num_Row
        temp_Signal		=	zeros( M, num_Col, 'like', X );
        temp_Signal( 1 : num_Row, : )	=	X;
    else
        temp_Signal	=	X( 1 : M, : );
        str_Warning_Msg	=	strcat( mfilename, ': The number of DAFT points is smaller than the number of rows of the input_signal' );
        warning( str_Warning_Msg );
    end


    chirp_Index	=	( 0 : M - 1 ).';
    Lamda1 = exp( 1i*2*pi*c1 * chirp_Index.^2 );
    Lamda2 = exp( 1i*2*pi*c2 * chirp_Index.^2 );

    X_temp1=repmat(Lamda2, 1, num_Col).*X;
    X_temp2=ifft(X_temp1).*sqrt(M);
    X_temp3=repmat(Lamda1, 1, num_Col).*X_temp2;
    x=X_temp3;

end

end
