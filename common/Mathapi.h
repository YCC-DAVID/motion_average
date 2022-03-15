#pragma once

void sort2(unsigned int n, double* ra,unsigned int* rb)
{
	if (n <= 1)
		return;

	unsigned int l,ir,i,j;
	double rra;
	unsigned int rrb;
	l = n / 2 + 1; 
	ir = n;
	do
	{
		if (l > 1)
		{
			l = l - 1;
			rra = ra[l - 1];
			rrb = rb[l - 1];
		}
		else
		{
			rra = ra[ir - 1];
			rrb = rb[ir - 1];
			ra[ir - 1] = ra[0];
			rb[ir - 1] = rb[0];
			ir = ir - 1;
			if (ir == 1)
			{
				ra[1 - 1] = rra;
				rb[1 - 1] = rrb;
				return;
			}
		}
		i = l;
		j = l + l;
		while (j <= ir)
		{
			if (j <= ir)
			{
				if (j < ir)
				{
					if (ra[j - 1] < ra[j])
					{
						j = j + 1;
					}
				}
				if (rra < ra[j - 1])
				{
					ra[i - 1] = ra[j - 1];
					rb[i - 1] = rb[j - 1];
					i = j;
					j = j + j;
				}
				else
				{
					j = ir + 1;
				}
			}
		}
		ra[i - 1] = rra;
		rb[i - 1] = rrb;
	}while(1);

}