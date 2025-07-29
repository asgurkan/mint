# Snakemake Kütüphanesi ve Popülasyon Genetiği Analizi

Bu belge, popülasyon genetiği alanında yaygın olarak uygulanan analizleri gerçekleştirmek üzere geliştirilmiş bir Snakemake pipeline’ının kullanımını detaylı bir şekilde açıklamaktadır. VCF (Variant Call Format) dosyasını girdi olarak kabul eden bu analiz zinciri, biyoinformatik analizlerin standart, modüler ve yeniden uygulanabilir bir biçimde yürütülmesine olanak sağlar.
Bu dokümantasyon, Coraman LAB bünyesinde düzenlenen "Snakemake Kütüphanesi Çalıştayı" kapsamında ekip üyelerine yönelik olarak hazırlanmıştır. Amaç, yalnızca mevcut pipeline’ın çalıştırılması değil, aynı zamanda bir Snakemake altyapısının kurulması, temel konfigürasyonların tanımlanması, girdi dosyalarının hazırlanması ve analiz çıktılarının tartışılmasıdır.

Pipeline, aşağıdaki temel analizleri kapsamaktadır:

* PCA (Principal Component Analysis)

* ADMIXTURE analizi

* FST (Fixation Index) hesaplaması

Snakemake’in sunduğu modüler yapı ve otomatik kaynak yönetimi sayesinde, bu analizler esnek biçimde yeniden çalıştırılabilir, farklı projelere kolayca uyarlanabilir ve çıktıların tutarlılığı güvence altına alınabilir.

## Sistem Gereksinimleri ve Ortam Kurulumu

Bu pipeline, Python ve Snakemake üzerine inşa edilmiştir. Öncelikle temel bir Conda ortamının oluşturulması gerekmektedir:

### Ortam kurulumu (tek bir temel ortam yeterlidir):

```bash
conda env create -f satup.yaml
conda activate mint_env
```

> Not: Snakemake kurulduğunda, kendine ait kurallar için gereken Conda ortamlarını tanımlı `envs/*.yaml` dosyaları üzerinden otomatik olarak oluşturacaktır.

Ardından pipeline girdisi olan ve "Örnek", "Popülasyon", "Proje" ve "Path" bilgilerini içeren sample.tsv dosyasının oluşturulması için hazırlanan bash script çalıştırılmalıdır.
```bash
chmod +x sample_file_generation.sh
./sample_file_generation.sh <vcf_file> <project_name>
```
---

## Proje Dosya Yapısı

```bash
project_directory/
├── config.yaml
├── sample.tsv
├── Snakefile
├── envs/
├── scripts/
├── data/
│   └── <project_name>/
└── logs/
```

---

## Girdi Dosyalarının Açıklamaları

### config.yaml

Pipeline’ın genel parametrelerini içerir.

**Örnek:**

```yaml
fst:
  window-size: 50000
  window-step: 25000
```

* `window-size`: FST analizinde kullanılacak pencere uzunluğu (örneğin: 50.000 bp)
* `window-step`: Pencereler arası kayma miktarı

Bu parametreler, FST analizinin gerçekleştirildiği kural içerisindeki parametre bölümünde 
```
params: 
    window_size = config["fst"]["window-size"]
```
 şeklinde çağırılır ve ilgili parametre shell komutu veya script içerisinde dinamik olarak çağırılır. 

### sample.tsv

Projeye ait örneklerin bilgisini barındırır. Minimum olarak proje ismi ve VCF dosyasının yolu verilmelidir.

**Örnek:**

```tsv
project	vcf_path
mysdav	data/mysdav/variants/filtered.vcf
```

| Sütun Adı  | Açıklama                                               |
| ---------- | ------------------------------------------------------ |
| `project`  | Proje adı; çıktı klasörleri bu isimle organize edilir. |
| `vcf_path` | Analiz edilecek .vcf dosyasının yolu                   |

---

## Kullanım

1. `config.yaml` ve `sample.tsv` dosyalarının projeye özel olarak düzenlenmesi.
2. Snakemake ortamının aktive edilmesi:

```bash
conda activate mint_env
```

3. Pipeline’ı çalıştırın:

```bash
snakemake --use-conda --conda-frontend conda -j 4 
```

---

## Analizler

### PCA (Principal Component Analysis)

PCA analizi, bireyler arası genetik varyasyonu özetlemek ve populasyon yapısını görselleştirmek için kullanılır. Bu adımda `plink` yazılımı aracılığıyla, LD (Linkage Disequilibrium) pruning işleminden geçmiş SNP verisi üzerinde özvektör (eigenvector) ve özdeğer (eigenvalue) hesaplaması gerçekleştirilir.

- **Girdi:** LD-pruned `.bed`, `.bim`, `.fam` dosyaları  
- **Kullanılan araç:** `plink --pca`  
- **Çıktı:**
  - `pca.eigenvec`: Örneklerin ana bileşenlere göre konumlarını içerir
  - `pca.eigenval`: Her bileşenin açıklanan varyans oranı  
- **Görselleştirme:** İlk iki bileşene göre bireylerin konumlandığı PCA scatter plot, örnek metadatalarıyla (örneğin popülasyon, grup vs.) renklendirilir

### ADMIXTURE Analizi

ADMIXTURE, bireylerin genetik bileşimlerini tahmin ederek her bir bireyin kaç farklı ata popülasyonundan geldiğini K bileşenli bir modelle hesaplar. Her `K` değeri için çapraz doğrulama hatası (CV error) raporlanır.

- **Girdi:** `.bed`, `.bim`, `.fam` dosyaları  
- **Kullanılan araç:** `admixture --cv <K>`  
- **Çıktı:**
  - `pruned.K.Q`: Bireylerin her ata popülasyona ait olasılıkları (Q-matrix)
  - `pruned.K.P`: Her ata popülasyonun SNP frekansları (P-matrix)
  - `admixture_cv_error_K{K}.txt`: K’ya karşılık gelen CV hatası  
- **Görselleştirme:** CV hata grafiği (K değerine göre hata) ve Q matrisinin barplot görselleştirmesi

### FST Analizi

FST, popülasyonlar arasındaki genetik farklılıkların derecesini ölçmek için kullanılır. Bu analizde `vcftools` ile pencere bazlı FST hesaplaması yapılır.

- **Girdi:** `.vcf` dosyası ve her popülasyona ait birey listesi (plain text)  
- **Kullanılan araç:** `vcftools --weir-fst-pop`  
- **Parametreler:**
  - `--fst-window-size`: Her pencerenin uzunluğu (örnek: 50.000 bp)
  - `--fst-window-step`: Pencereler arası kayma miktarı  
- **Çıktı:** `.windowed.weir.fst` dosyası  
- **Görselleştirme:** Manhattan tipi FST plotu ile genomik konumlara göre FST dağılımı

### Raporlama

Tüm analizlerden elde edilen çıktılar tek bir özetleme adımında toplanır. Bu adım, downstream analizleri kolaylaştırmak, raporlama sürecini otomatize etmek ve kalite kontrol sürecine katkı sunmak için uygulanır.

- **Girdi:**
  - PCA görselleri ve dosyaları
  - ADMIXTURE çıktıları
  - FST analiz çıktısı ve görselleştirmesi
  - Annotation dağılımları ve ek analiz klasörleri  
- **Çıktı:**
  - `.txt` formatında sadeleştirilmiş, yorumlanabilir proje özeti
  - Dosya yolları, kullanılan parametreler ve görsellerin konumları dahil edilir

---

## Örnek Çıktı Yapısı

```
data/
├── mysdav/
│   ├── pca/
│   │   ├── pca.eigenvec
│   │   └── pca.eigenval
│   ├── admixture/
│   │   ├── K3/
│   │   │   ├── pruned.3.Q
│   │   │   ├── pruned.3.P
│   │   ├── K4/
│   │   │   ├── pruned.4.Q
│   │   │   ├── pruned.4.P
│   │   ├── K5/
│   │   │   ├── pruned.5.Q
│   │   │   ├── pruned.5.P
│   ├── fst/
│   │   ├── mysdav.windowed.weir.fst
│   └── figures/
│       ├── pca/pca_plot.png
│       ├── admixture/error_plot.png
│       └── fst/fst_manhattan_locus.png
```

---

## Hata Giderme

* Hatalarla karşılaşılması durumunda ilk olarak `logs/` klasöründeki log dosyalarının kontrol edilmesi
* Anlamlandırılmayan hatalar için kurulan snakemake environment dosyasının aktive edilerek ilgili rule'un başka bir terminalde çalıştırılarak çıktılarının kontrol edilmesi
* Kullanılan kütüphanelerin versiyonlarının güncel ve uyumlu olanlarla değiştirilmesi

önerilmektedir. 

---
Soru ve önerileriniz için :
ahmetsametgurkan@gmail.com

İyi çalışmalar.
