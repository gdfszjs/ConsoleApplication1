#pragma once

#include "tuple.h"
#include "Array.h"
#include "PBCG.h"

using namespace LibIV::Math;
using namespace LibIV::Memory::Array;

namespace LibIV
{
	namespace Poisson
	{
		struct INFO
		{
			int N;    // The number of neighbors
			v3i self; // (x,y,index)
			v4i left; // (isExist, x, y, index)
			v4i right;
			v4i up;
			v4i down;
		};

		inline double PIE_b(const INFO & pix, const TensorInt & src, const TensorInt & dst)
		{
			double ans = 0;

			int width = src.cols();
			int height = src.rows();

			ans += pix.N * src(pix.self[0],pix.self[1]);

			/*if(pix.left[0])
			{*/
			ans -= src(pix.left[1],pix.left[2]);

			if(pix.left[1] == 0)
				ans += src(pix.self[0],pix.self[1]);
			else if(pix.left[3] == -1)
				ans += dst(pix.left[1],pix.left[2]);
			//}

			/*if(pix.right[0])
			{*/
			ans -= src(pix.right[1],pix.right[2]);

			if(pix.right[1] == width - 1)
				ans += src(pix.self[0],pix.self[1]);
			else if(pix.right[3] == -1)
				ans += dst(pix.right[1],pix.right[2]);
			//}

			/*if(pix.up[0])
			{*/
			ans -= src(pix.up[1],pix.up[2]);

			if(pix.up[2] == 0)
				ans += src(pix.self[0],pix.self[1]);
			else if(pix.up[3] == -1)
				ans += dst(pix.up[1],pix.up[2]);
			//}

			/*if(pix.down[0])
			{*/
			ans -= src(pix.down[1],pix.down[2]);

			if(pix.down[2] == height-1)
				ans += src(pix.self[0],pix.self[1]);
			else if(pix.down[3] == -1)
				ans += dst(pix.down[1],pix.down[2]);
			//}

			return ans;
		}

		inline void possionImageEditing_one_channel(TensorByte & ans, const TensorInt & src, const TensorInt & dst, const TensorByte & msk)
		{
			int width = src.cols(), height = src.rows();

			_check_(dst.rows() == height && dst.cols() == width && msk.rows() == height && msk.cols() == width);

			ans.fill(0);

			TensorInt index;
			index.set(width,height);
			index.fill(-1);

			// Get index
			int pixel_num = 0;
			for(int i = 0;i<height;i++)
			{
				for(int j = 0;j<width;j++)
				{
					if(msk(j,i))
						index(j,i) = pixel_num++;
				}

			}

			if(pixel_num == 0) return;

			// parsing the connections between pairs of pixels
			INFO * vars = new INFO[pixel_num];

			pixel_num = 0;

			for(int i = 0;i<height;i++)
			{
				for(int j = 0;j<width;j++)
				{
					if(msk(j,i))
					{
						vars[pixel_num].N = 0;
						vars[pixel_num].left[0] = 0;
						vars[pixel_num].right[0] = 0;
						vars[pixel_num].up[0] = 0;
						vars[pixel_num].down[0] = 0;


						// self
						vars[pixel_num].self[0] = j;
						vars[pixel_num].self[1] = i;
						vars[pixel_num].self[2] = index(j,i);

						// neighbors
						/*if(j > 0)
						{*/
						vars[pixel_num].N++;
						vars[pixel_num].left[0] = 1;
						vars[pixel_num].left[1] = j-1;
						vars[pixel_num].left[2] = i;
						vars[pixel_num].left[3] = index(j-1,i);
						//}

						/*if(j < width - 1)
						{*/
						vars[pixel_num].N++;
						vars[pixel_num].right[0] = 1;
						vars[pixel_num].right[1] = j+1;
						vars[pixel_num].right[2] = i;
						vars[pixel_num].right[3] = index(j+1,i);
						//}

						/*if(i > 0)
						{*/
						vars[pixel_num].N++;
						vars[pixel_num].up[0] = 1;
						vars[pixel_num].up[1] = j;
						vars[pixel_num].up[2] = i-1;
						vars[pixel_num].up[3] = index(j,i-1);

						//}

						/*if(i < height - 1)
						{*/
						vars[pixel_num].N++;
						vars[pixel_num].down[0] = 1;
						vars[pixel_num].down[1] = j;
						vars[pixel_num].down[2] = i+1;
						vars[pixel_num].down[3] = index(j,i+1);
						//}
						pixel_num++;
					}
				}
			}

			// construct sparse matrix
			double * sa = new double[5 * pixel_num];
			unsigned long * ija = new unsigned long[5 * pixel_num];
			double * x = new double[pixel_num];
			double * b = new double[pixel_num];

			memset(sa,0,sizeof(double) * pixel_num * 5);
			memset(ija,0,sizeof(unsigned long) * pixel_num * 5);
			memset(x,0,sizeof(double)*pixel_num);
			memset(b,0,sizeof(double)*pixel_num);

			int _i = pixel_num + 1;
			for(int i = 0;i<pixel_num;i++)
			{
				sa[i] = vars[i].N;
				ija[i] = _i;

				if(/*vars[i].up[0] && */vars[i].up[3] != -1)
				{
					ija[_i] = vars[i].up[3];
					sa[_i++] = -1;
				}

				if(/*vars[i].left[0] && */vars[i].left[3] != -1)
				{
					ija[_i] = vars[i].left[3];
					sa[_i++] = -1;
				}

				if(/*vars[i].right[0] && */vars[i].right[3] != -1)
				{
					ija[_i] = vars[i].right[3];
					sa[_i++] = -1;
				}

				if(/*vars[i].down[0] && */vars[i].down[3] != -1)
				{
					ija[_i] = vars[i].down[3];
					sa[_i++] = -1;
				}

				b[i] = PIE_b(vars[i],src,dst);
			}
			ija[pixel_num] = _i;


			int mpcgx_iter;
			double mpcgx_error;
			PBCG mpcgx(sa,ija);
			mpcgx.linbcg(x,b,pixel_num,mpcgx_iter,mpcgx_error,3,1e-5,200);

			pixel_num = 0;
			for(int i = 0;i<height;i++)
			{
				for(int j = 0;j<width;j++)
					if(msk(j,i))
					{
						if(x[pixel_num] < 0)
							ans(j,i) = 0;
						else if(x[pixel_num] > 255)
							ans(j,i) = (uchar)255;
						else
							ans(j,i) = (uchar)(x[pixel_num]);
						pixel_num++;
					}
			}

			delete [] sa;
			delete [] ija;
			delete [] x;
			delete [] b;
			delete [] vars;
		}

		// Attention:
		// dst_img->nChannels can be 3 or 4
		// src_pixel: the non-existing pixel should be filled with (-1,-1,-1)
		inline void possionImageEditing(TensorV3i & src_pixel, IplImage * dst_img)
		{
			if(!dst_img) return;

			int width = dst_img->width ,height = dst_img->height;

			_check_(src_pixel.rows() == height && src_pixel.cols() == width);

			// Get MSK
			TensorByte msk;
			msk.set(width,height);
			msk.fill(0);

			for(int i = 0;i<height;i++)
			{
				for(int j = 0;j<width;j++)
				{
					if(src_pixel(j,i)[0] != -1)
					{
						bool isEdge = false;

						if(j == 0 || src_pixel(j-1,i)[0] == -1)
							isEdge = true;

						if(j == width-1 || src_pixel(j+1,i)[0] == -1)
							isEdge = true;

						if(i == 0 || src_pixel(j,i-1)[0] == -1)
							isEdge = true;

						if(i == height-1 || src_pixel(j,i+1)[0] == -1)
							isEdge = true;

						if(!isEdge)
							msk(j,i) = 1;
					}
				}
			}

			TensorByte ans_r(height,width),ans_g(height,width),ans_b(height,width);
			TensorInt src_r(height,width),src_g(height,width),src_b(height,width);
			TensorInt dst_r(height,width),dst_g(height,width),dst_b(height,width);
			int r,g,b;

			for(int i = 0;i<height;i++)
			{
				for(int j = 0;j<width;j++)
				{	
					src_b(j,i) = src_pixel(j,i)[0];
					src_g(j,i) = src_pixel(j,i)[1];
					src_r(j,i) = src_pixel(j,i)[2];

					LibIV::CppHelpers::Global::imgData(r,g,b,dst_img,j,i);

					dst_r(j,i) = r;
					dst_g(j,i) = g;
					dst_b(j,i) = b;
				}
			}

			possionImageEditing_one_channel(ans_r,src_r,dst_r,msk);
			possionImageEditing_one_channel(ans_g,src_g,dst_g,msk);
			possionImageEditing_one_channel(ans_b,src_b,dst_b,msk);

			for(int i = 0;i<height;i++)
			{
				for(int j = 0;j<width;j++)
					if(msk(j,i))
						LibIV::CppHelpers::Global::setImgData(ans_r(j,i),ans_g(j,i),ans_b(j,i),dst_img,j,i);
			}
		}

	}
}