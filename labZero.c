
#include "utils.h"
#include "ZeroFuncao.h"

int main() {
    fesetround(FE_DOWNWARD); // Define arredondamento para baixo
    Polinomio pol;
    real_t raiz, a, b, x0;
    int it = 0;
    rtime_t tempo; 

    scanf("%d", &pol.grau); // Grau do polinomio

    pol.p = (real_t*) malloc(sizeof(real_t) * (pol.grau + 1));

    if (pol.p == NULL) {
        printf("Erro ao alocar memoria para os coeficientes do polinomio\n");
        return 1;
    }

    for (int i = pol.grau; i >= 0; --i) scanf("%lf", &pol.p[i]); // Coeficientes do polinomio

    scanf("%lf %lf", &a, &b); // intervalo onde esta uma das raizes

    x0 = (a + b) / 2;

    printf("\nRAPIDO\n\n");
    
    for(int i = 1; i < 4; i++){
        tempo = timestamp();
        bisseccao(pol, a, b, i, &it, &raiz, &calcPolinomio_rapido);
        tempo = timestamp() - tempo;
        printf("bissec %.15e valor_crit-0%i %i %.8e\n", raiz, i, it, tempo);
    }

    for(int i = 1; i < 4; i++){
        tempo = timestamp();
        newtonRaphson(pol, x0, i, &it, &raiz, &calcPolinomio_rapido);
        tempo = timestamp() - tempo;
        printf("newton %.15e valor_crit-0%i %i %.8e\n", raiz, i, it, tempo);
    }

    
    printf("\nLENTO\n\n");

    for(int i = 1; i < 4; i++){
        tempo = timestamp();
        bisseccao(pol, a, b, i, &it, &raiz, &calcPolinomio_lento);
        tempo = timestamp() - tempo;
        printf("bissec %.15e valor_crit-0%i %i %.8e\n", raiz, i, it, tempo);
    }

    for(int i = 1; i < 4; i++){
        tempo = timestamp();
        newtonRaphson(pol, x0, i, &it, &raiz, &calcPolinomio_lento);
        tempo = timestamp() - tempo;
        printf("newton %.15e valor_crit-0%i %i %.8e\n", raiz, i, it, tempo);
    }

    free(pol.p);

    return 0;
}
