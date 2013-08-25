Bayesovská regrese
===

Potřebné python moduly
---
- numpy
- qt4

Popis
---

Program demonstruje bayesovskou regresi na jednoduchém příkladu, kde je možno nastavit několik parametrů a sledovat výsledek na grafu.

Program nepracuje s reálnými daty, ale generuje náhodná data (při velikosti 30 by skutečná data stejně byla podobná náhodným).

Rozhodl jsem se, že všechny "módy" budou mít 50 bázových funkcí. To se neukázalo jako dobrá volba.

Přiznám se totiž, že nemám úplně vhled do toho, jak kovarianční matice konkrétně ovlivňuje které proměnné. Pokud se potom např. dají bázové funkce jako mocniny od x^0 do x^49, je nutné být u velkých exponentů velmi opatrný - x^49 roste výrazně rychleji, než např x^5, a už při malém kladném koeficientu začne přerůstat všechno ostatní.

Kromě toho velké množství funkcí má tendenci overfitovat.

V grafu jsou vždy vykresleny trénovací data - pro x od 1 do 30 je zde trénovací příklad, znázorněný modrým kolečkem. Jak bude popsáno dále, data jsou vybírána náhodně, přičemž při každé nové volbě se vytáhnou úplně nová. To sice mírně zhoršuje srovnatelnost, ale na druhou stranu je vidět "více experimentů".

Možnosti
---

### Možnosti zobrazení

Je možno zobrazit bázové funkce, náhodné priory a náhodné posteriory.

Bázové funkce nemá např. u sinů tolik cenu kreslit, protože se nakreslí všechny naráz. 

Náhodné priory a náhodné posteriory se samplují z normálních a normal-gamma rozdělení. Tyto se počítají podle už odevzdaných vzorců. 

### Kovariance vah na prioru

Jak jsem uvedl výše, schází mi trochu "vhled" do kovariančních matic a jejich významu.

Je možné vybrat z několika možností:

- jednotková matice
- "kolem jednotkové" matice - na diagonále jsou největší hodnoty, ve vzálenosti 1 od diagonály o něco menší, ve vzálenosti 2 jedničky, všude jinde nuly
- "vlevo nahoře" - v matici 17x17 vlevo nahoře jsou samé jedničky (krom diagonály, kde jsou dvojky), zbytek nuly (krom diagonály, kde je 10e-30)
- "vpravo dole" - totéž, jen jsou jedničky vpravo dole
- "rnd" - náhodná kovarianční matice, hodnoty 0-1
- "x10" - vynásobí kov. matici deseti
- "x100" - vynásoví stem (oboje zaškrtnuté dohromady vynásobí tisícem)

### Střední hodnota vah prioru

Prioru vah se můžou dát různé hodnoty.

- vše 0
- vše 1000
- začátek (prvních 17) tisíc, zbytek 0
- totéž s koncem

### Funkce

Mám tři množiny bázových funkcí.

- polynomy - sečítám x^(n/10)
- sinus posun - sečítám sin(x+(2*pi/50)). Tj. siny se stejnou frekvencí, jenom různě posunuté.
- sinus posun+frekvence - sečítám siny jednak posunuté, ale také s 3 různými frekvencemi, přičemž váhy na začátku jsou k nižším frekvencím a váhy na konci k vyšším

U sinů mám navíc to, že pro n==0 vracím konstantní funkci 1.

### Data

Jak generuji trénovací data.

- náhodné kolem nuly
- náhodné kolem 30 (se stejnou odchylkou)
- absolutní hodnota z náhodné kolem 1 * sinus [tím chci simulovat to, že správná distribuce je "skoro sinus", ačkoliv není přesně podle modelu]
- totéž, jen je sinus 2x rychlejší

### Známá odchylka

Dále je možné zvolit známou odchylku v předpokládaném modelu, nebo parametr alfa pro gamma distribuci odchylky. (Beta parametr je konstantní.)

Zvolil jsem rozdílné hodnoty, aby byl víc vidět rozdíl.
