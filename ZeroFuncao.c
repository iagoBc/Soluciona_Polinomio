
#include "ZeroFuncao.h"


int ULPs(double A, double B) {
    Double_t dA, dB;

    // Converter os doubles para a estrutura Double_t
    memcpy(&dA.f, &A, sizeof(double));
    memcpy(&dB.f, &B, sizeof(double));

    // Sinais diferentes significam que eles não correspondem
    if (dA.parts.sign != dB.parts.sign) {
        return (dA.f == dB.f); // Verifique a igualdade para ter certeza de que +0==-0
    }

    // Calcula a diferença em ULPs
    int64_t dif_ulps = llabs(dA.i - dB.i);

    if (dif_ulps <= ULPS) return 1;
    return 0;
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, void f(Polinomio p, real_t x, real_t *fx, real_t *dfx))
{
    real_t x_novo = x0, x_antigo;
    real_t fx, dfx;
    (*it) = 0;

    do {
        // Calcula f(x) e f'(x) no ponto x_novo
        f(p, x_novo, &fx, &dfx);

        // Verifica se a derivada = zero (evitar divisao por zero)
        if (fabs(dfx) <= EPS) return x_novo;  // Retorna o valor atual de x_novo
        
        // Armazena o valor anterior de x_novo
        x_antigo = x_novo;

        // Atualiza o valor de x_novo
        x_novo = x_antigo - (fx / dfx);
        
        (*it)++;

        // Criterio-01: |x_novo − x_antigo| ⩽ EPS
        if (criterioParada == 1 && fabs(x_novo - x_antigo) <= EPS) {
            *raiz = x_novo;
            return x_novo;
        }

        // Criterio-02:  |f(x_novo)| ⩽ DBL_EPSILON
        if (criterioParada == 2 && fabs(fx) <= EPS) {
            *raiz = x_novo;
            return x_novo;
        }
        
        // Criterio-03: ULP's entre x_novo e  x_antigo ⩽ 2
        if (criterioParada == 3 && ULPs(x_novo, x_antigo)) {
            *raiz = x_novo;
            return x_novo;
        }

    } while (*it <= MAXIT);  // Limite de iteracoes

    *raiz = x_novo;
    return x_novo;
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, void f(Polinomio p, real_t x, real_t *fx, real_t *dfx))
{
    double x_novo, x_antigo;
    real_t fa, fx, dpx; // fa = f(a), fx = f(x_novo)

    x_novo = (a + b)/2;

    (*it) = 0;

    f(p, x_novo, &fx, &dpx);
    f(p, a, &fa, &dpx);

    // Verifica se a raiz está em [a, x_novo] ou [x_novo, b]
    if(fa * fx < 0) b = x_novo;  // Se f(a) * f(x_novo) < 0, a raiz está entre a e x_novo, então ajusta b
    else if(fa * fx > 0) a = x_novo;  // Se f(a) * f(x_novo) > 0, a raiz está entre x_novo e b, então ajusta a
    else{
        *raiz = x_novo;  // Se f(a) * f(x_novo) == 0, então x_novo é a raiz
        return x_novo;
    }

    do{    
        // Salva a posição anterior de x_novo
        x_antigo = x_novo;
        
        // Atualiza x_novo para o novo ponto médio
        x_novo = (a + b) / 2;
        
        // Recalcula os valores de f(x_novo) e f(a)
        f(p, x_novo, &fx, &dpx);
        f(p, a, &fa, &dpx);

        (*it)++;

        // Verifica novamente em qual subintervalo a raiz se encontra
        if(fa * fx < 0) b = x_novo;  // Ajusta b para o novo intervalo
        else if(fa * fx > 0) a = x_novo;  // Ajusta a para o novo intervalo
        else{
            *raiz = x_novo;  // Se f(a) * f(x_novo) == 0, então x_novo é a raiz
            return x_novo;
        }

        // Criterio-01: |x_novo − x_antigo| ⩽ EPS
        if (criterioParada == 1 && fabs(x_novo - x_antigo) <= EPS) {
            *raiz = x_novo;
            return x_novo;
        }

        // Criterio-02:  |f(x_novo)| ⩽ DBL_EPSILON
        if (criterioParada == 2 && fabs(fx) <= EPS) {
            *raiz = x_novo;
            return x_novo;
        }
        
        // Criterio-03: ULP's entre x_novo e  x_antigo ⩽ 2
        if (criterioParada == 3 && ULPs(x_novo, x_antigo)) {
            *raiz = x_novo;
            return x_novo;
        }

    } while (*it <= MAXIT);  // Limite de iteracoes

    *raiz = x_novo;
    return x_novo;
}

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    double b = 0, c = 0;

    for(int i = p.grau; i > 0; --i){      
        b = b*x + p.p[i];   // Calcula o valor do polinômio f(x) de forma eficiente (avaliação Horner)
        c = c*x + b;        // Calcula a derivada f'(x) de forma eficiente (usando a fórmula de Horner)
    }

    
    b = b*x + p.p[0];
    *px = b;  // Armazena o valor do polinômio f(x)
    *dpx = c; // Armazena o valor da derivada f'(x)
}

void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    *px = 0;  
    *dpx = 0; 

    for(int i = 0; i <= p.grau; i++){
        *px += p.p[i] * pow(x,i); // Calcula o valor do polinômio de forma direta (potência de x)
        if (i > 0) {
            *dpx += i * p.p[i] * pow(x, i - 1);  // Derivada usando a regra do expoente
        }
    }
}

