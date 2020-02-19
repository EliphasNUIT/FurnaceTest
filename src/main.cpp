#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "dds.h"

float const MATH_PI = 3.14159f;

unsigned const sampleNum = 8192;

void ExportTable(char const* name, float* lutData, unsigned lutWidth, unsigned lutHeight)
{
	FILE* f = fopen(name, "w");
	for (unsigned y = 0; y < lutHeight; ++y)
	{
		for (unsigned x = 0; x < lutWidth; ++x)
		{
			float const ndotv = (x + 0.5f) / lutWidth;
			float const roughness = (y + 0.5f) / lutHeight;
			fprintf(f, "%f, %f, %f\n", ndotv, roughness, lutData[x * 4 + y * lutWidth * 4]);
		}
	}
	fclose(f);
}

float Saturate(float x)
{
	return std::min(std::max(x, 0.0f), 1.0f);
}

uint32_t ReverseBits(uint32_t v)
{
	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
	v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
	v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
	v = (v >> 16) | (v << 16);
	return v;
}

struct HalfVector {
	float x;
	float y;
	float z;
};

HalfVector ImportanceGGX(float e1, float e2, float roughness) {
	struct HalfVector H;

	float a = roughness * roughness;
	float a2 = a * a;

	float Phi = 2.0f * MATH_PI * e1;
	float CosTheta = sqrtf((1.0f - e2) / (1.0f + (a2 - 1.0f) * e2));
	float SinTheta = sqrtf(1.0f - CosTheta * CosTheta);

	H.x = SinTheta * cosf(Phi);
	H.y = SinTheta * sinf(Phi);
	H.z = CosTheta;

	return H;
}

HalfVector ImportanceSampleEstevezKulla(float e1, float e2, float roughness) {
	struct HalfVector H;

	float a = roughness;

	float Phi = 2.0f * MATH_PI * e1;
	float SinTheta = powf(e2, a / (2.0f * a + 1.0f));
	float CosTheta = sqrtf(1.0f - SinTheta * SinTheta);

	H.x = SinTheta * cosf(Phi);
	H.y = SinTheta * sinf(Phi);
	H.z = CosTheta;

	return H;
}

float vMax(float a, float b) {
	return a >= b ? a : b;
}

float GeomGGX(float roughness, float NoV, float NoL) {
	float a = roughness * roughness;
	float a2 = a * a;
	float NoV2 = vMax(NoV * NoV, 1e-6f);
	float NoL2 = vMax(NoL * NoL, 1e-6f);
	float Vis_SmithV = sqrtf(1.0f - a2 + a2 / NoV2);
	float Vis_SmithL = sqrtf(1.0f - a2 + a2 / NoL2);
	return 2.0f / (Vis_SmithV + Vis_SmithL);
};

float NoLVisPDFGGX(float roughness, float NoV, float NoL, float VoH, float NoH) {
	float GVis = 0.25f * GeomGGX(roughness, NoV, NoL) / vMax(NoL * NoV, 1e-6f);
	float NoLPDF = 4.0f * NoL * VoH / vMax(NoH, 1e-6f);
	return NoLPDF * GVis;
};

float NoLVisPDFAshikminG(float roughness, float NoV, float NoL, float VoH, float NoH) {
	float GVis = 0.25f / vMax(NoL + NoV - NoV * NoL, 1e-6f);
	float NoLPDF = 4.0f * NoL * VoH / vMax(NoH, 1e-6f);
	return NoLPDF * GVis;
};

float L(float cosTheta, float roughness) {
	float mixCst = (1.0f - roughness) * (1.0f - roughness);
	double a = 25.3245 * mixCst + 21.5473 * (1.0 - mixCst);
	double b = 3.32435 * mixCst + 3.82987 * (1.0 - mixCst);
	double c = 0.16801 * mixCst + 0.19823 * (1.0 - mixCst);
	double d = -1.27393 * mixCst + -1.97760 * (1.0 - mixCst);
	double e = -4.85967 * mixCst + -4.32054 * (1.0 - mixCst);
	return (float)(a / (1.0 + b * pow(cosTheta, c)) + d * cosTheta + e);
}

float Lambda(float cosTheta, float roughness) {
	if (abs(cosTheta) < 0.5f) {
		return expf(L(abs(cosTheta), roughness));
	}
	return expf(2.0f * L(0.5f, roughness) - L(1.0f - abs(cosTheta), roughness));
}

float EstevezKullaG(float roughness, float NoV, float NoL) {

	return 1.0f / vMax(1.0f + Lambda(NoV, roughness) + Lambda(NoL, roughness), 1e-6f);
}

float NoLVisPDFEstevezKullaG(float roughness, float NoV, float NoL, float VoH, float NoH) {
	float GVis = 0.25f * EstevezKullaG(roughness, NoV, NoL) / vMax(NoV * NoL, 1e-6f);
	float NoLPDF = 4.0f * NoL * VoH / vMax(NoH, 1e-6f);
	return NoLPDF * GVis;
};

float WhiteFurnaceTest_GGX(float roughness, float ndotv)
{

	float const vx = sqrtf(1.0f - ndotv * ndotv);
	float const vy = 0.0f;
	float const vz = ndotv;

	float integral = 0.0f;
	struct HalfVector H;
	for (unsigned i = 0; i < sampleNum; ++i)
	{
		float const e1 = (float)i / sampleNum;
		float const e2 = (float)((double)ReverseBits(i) / (double)0x100000000LL);

		H = ImportanceGGX(e1, e2, roughness);

		float const hx = H.x;
		float const hy = H.y;
		float const hz = H.z;

		float const vdothUnsat = vx * hx + vy * hy + vz * hz;
		float const lx = 2.0f * vdothUnsat * hx - vx;
		float const ly = 2.0f * vdothUnsat * hy - vy;
		float const lz = 2.0f * vdothUnsat * hz - vz;

		float const ndotl = std::max(lz, 0.0f);
		float const ndoth = std::max(hz, 0.0f);
		float const vdoth = std::max(vdothUnsat, 0.0f);

		integral += NoLVisPDFGGX(roughness, ndotv, ndotl, vdoth, ndoth);
	}
	integral *= 0.5 / sampleNum;
	return integral;
}

float WhiteFurnaceTest_EK_AG(float roughness, float ndotv)
{

	float const vx = sqrtf(1.0f - ndotv * ndotv);
	float const vy = 0.0f;
	float const vz = ndotv;

	float integral = 0.0f;
	struct HalfVector H;
	for (unsigned i = 0; i < sampleNum; ++i)
	{
		float const e1 = (float)i / sampleNum;
		float const e2 = (float)((double)ReverseBits(i) / (double)0x100000000LL);

		H = ImportanceSampleEstevezKulla(e1, e2, roughness);

		float const hx = H.x;
		float const hy = H.y;
		float const hz = H.z;

		float const vdothUnsat = vx * hx + vy * hy + vz * hz;
		float const lx = 2.0f * vdothUnsat * hx - vx;
		float const ly = 2.0f * vdothUnsat * hy - vy;
		float const lz = 2.0f * vdothUnsat * hz - vz;

		float const ndotl = std::max(lz, 0.0f);
		float const ndoth = std::max(hz, 0.0f);
		float const vdoth = std::max(vdothUnsat, 0.0f);

		integral += NoLVisPDFAshikminG(roughness, ndotv, ndotl, vdoth, ndoth);
	}
	integral *= 0.5 / sampleNum;
	return integral;
}

float WhiteFurnaceTest_EK(float roughness, float ndotv)
{

	float const vx = sqrtf(1.0f - ndotv * ndotv);
	float const vy = 0.0f;
	float const vz = ndotv;

	float integral = 0.0f;
	struct HalfVector H;
	for (unsigned i = 0; i < sampleNum; ++i)
	{
		float const e1 = (float)i / sampleNum;
		float const e2 = (float)((double)ReverseBits(i) / (double)0x100000000LL);

		H = ImportanceSampleEstevezKulla(e1, e2, roughness);

		float const hx = H.x;
		float const hy = H.y;
		float const hz = H.z;

		float const vdothUnsat = vx * hx + vy * hy + vz * hz;
		float const lx = 2.0f * vdothUnsat * hx - vx;
		float const ly = 2.0f * vdothUnsat * hy - vy;
		float const lz = 2.0f * vdothUnsat * hz - vz;

		float const ndotl = std::max(lz, 0.0f);
		float const ndoth = std::max(hz, 0.0f);
		float const vdoth = std::max(vdothUnsat, 0.0f);

		integral += NoLVisPDFEstevezKullaG(roughness, ndotv, ndotl, vdoth, ndoth);
	}
	integral *= 0.5/sampleNum;
	return integral;
}

int main()
{
	unsigned const LUT_WIDTH = 128;
	unsigned const LUT_HEIGHT = 128;

	float lutWhiteFurnace_GGX[LUT_WIDTH * LUT_HEIGHT * 4];
	float lutWhiteFurnace_EK_AG[LUT_WIDTH * LUT_HEIGHT * 4];
	float lutWhiteFurnace_EK[LUT_WIDTH * LUT_HEIGHT * 4];

	for (unsigned y = 0; y < LUT_HEIGHT; ++y)
	{
		float const roughness = (y + 0.5f) / LUT_HEIGHT;

		for (unsigned x = 0; x < LUT_WIDTH; ++x)
		{
			float const ndotv = (x + 0.5f) / LUT_WIDTH;

			float const whiteFurnace_GGX = WhiteFurnaceTest_GGX(roughness, ndotv);

			lutWhiteFurnace_GGX[x * 4 + y * LUT_WIDTH * 4 + 0] = whiteFurnace_GGX;
			lutWhiteFurnace_GGX[x * 4 + y * LUT_WIDTH * 4 + 1] = whiteFurnace_GGX;
			lutWhiteFurnace_GGX[x * 4 + y * LUT_WIDTH * 4 + 2] = whiteFurnace_GGX;
			lutWhiteFurnace_GGX[x * 4 + y * LUT_WIDTH * 4 + 3] = 1.0f;


			float const whiteFurnace_EK_AG = WhiteFurnaceTest_EK_AG(roughness, ndotv);

			lutWhiteFurnace_EK_AG[x * 4 + y * LUT_WIDTH * 4 + 0] = whiteFurnace_EK_AG;
			lutWhiteFurnace_EK_AG[x * 4 + y * LUT_WIDTH * 4 + 1] = whiteFurnace_EK_AG;
			lutWhiteFurnace_EK_AG[x * 4 + y * LUT_WIDTH * 4 + 2] = whiteFurnace_EK_AG;
			lutWhiteFurnace_EK_AG[x * 4 + y * LUT_WIDTH * 4 + 3] = 1.0f;


			float const whiteFurnace_EK = WhiteFurnaceTest_EK(roughness, ndotv);

			lutWhiteFurnace_EK[x * 4 + y * LUT_WIDTH * 4 + 0] = whiteFurnace_EK;
			lutWhiteFurnace_EK[x * 4 + y * LUT_WIDTH * 4 + 1] = whiteFurnace_EK;
			lutWhiteFurnace_EK[x * 4 + y * LUT_WIDTH * 4 + 2] = whiteFurnace_EK;
			lutWhiteFurnace_EK[x * 4 + y * LUT_WIDTH * 4 + 3] = 1.0f;

			printf(".");
		}
	}

	ExportTable("assets/whiteFurnace_GGX.txt", lutWhiteFurnace_GGX, LUT_WIDTH, LUT_HEIGHT);
	ExportTable("assets/whiteFurnace_EK_AG.txt", lutWhiteFurnace_EK_AG, LUT_WIDTH, LUT_HEIGHT);
	ExportTable("assets/whiteFurnace_EK.txt", lutWhiteFurnace_EK, LUT_WIDTH, LUT_HEIGHT);

	SaveDDS("assets/whiteFurnace_GGX.dds", DDS_FORMAT_R32G32B32A32_FLOAT, 16, LUT_WIDTH, LUT_HEIGHT, lutWhiteFurnace_GGX);
	SaveDDS("assets/whiteFurnace_EK_AG.dds", DDS_FORMAT_R32G32B32A32_FLOAT, 16, LUT_WIDTH, LUT_HEIGHT, lutWhiteFurnace_EK_AG);
	SaveDDS("assets/whiteFurnace_EK.dds", DDS_FORMAT_R32G32B32A32_FLOAT, 16, LUT_WIDTH, LUT_HEIGHT, lutWhiteFurnace_EK);
	return 0;
}