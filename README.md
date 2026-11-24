# Network Reliability – Course Project

פרויקט גמר בקורס **אמינות רשתות**, תואר שני בהנדסת תוכנה (SCE).

**מגישים:** Kobi Loya, Roman Vikhman

## חלק א' – הערכת אמינות בשיטת Crude Monte Carlo

הערכת האמינות של רשת הדודקהדרון (20 קודקודים אמינים, 30 צלעות לא אמינות) בקריטריון
Terminal Connectivity, עם הטרמינלים **T1=4, T2=10, T3=17**.

- מצב הרשת בכל ניסוי נבדק באמצעות DSS (Union-Find) עם union by rank ו-path compression.
- הפתרון המוגש (`cmc_reliability.py`) כולל CMC בלבד, כנדרש בעבודה.
- בנפרד, ולא כחלק מההגשה, `validation/exact_validation.py` מחשב את האמינות **המדויקת** R(p) (frontier-based dynamic programming) כדי לבדוק את אומדני ה-CMC.

### מבנה התיקייה

```
part-a-cmc/
├── assignment/
│   └── project_part_a.pdf  # דף העבודה של חלק א'
├── cmc_reliability.py      # הפתרון המוגש: CMC + DSS
├── validation/
│   └── exact_validation.py # חישוב מדויק לאימות (לא חלק מההגשה)
├── results/
│   ├── run_output.txt                # פלט הפתרון (seed=2026)
│   └── exact_validation_output.txt   # פלט החישוב המדויק
└── report/
    ├── אמינות_רשתות_פרויקט_חלק_א.pdf    # הדוח להגשה
    ├── אמינות_רשתות_פרויקט_חלק_א.docx   # הדוח בפורמט Word
    ├── report.md           # מקור הדוח (Markdown + LaTeX)
    ├── network.png         # שרטוט הרשת
    ├── build.py            # בניית ה-PDF וה-DOCX מהמקור
    └── reference.docx      # תבנית עיצוב ל-Word
```

### הרצה

נדרש Python 3.8 ומעלה, ללא ספריות חיצוניות.

```bash
cd part-a-cmc
python cmc_reliability.py
```

ההרצה המלאה (כל ערכי p עם M=1000 ו-M=10000, וחמש ההרצות של סעיף E) לוקחת כשנייה.
להרצת האימות: `python validation/exact_validation.py`.
כדי לשנות טרמינלים, יש לערוך את השורה `TERMINALS` בקובץ.

### בניית הדוח מחדש

נדרשים `pandoc` ו-Google Chrome:

```bash
cd part-a-cmc/report
python build.py
```

### תוצאות

| p | M=1000 | M=10000 | R מדויק (אימות) |
|:---:|:---:|:---:|:---:|
| 0.1 | 0.0000 | 0.0000 | 0.000003 |
| 0.2 | 0.0000 | 0.0001 | 0.000370 |
| 0.3 | 0.0030 | 0.0069 | 0.006670 |
| 0.4 | 0.0610 | 0.0478 | 0.048590 |
| 0.5 | 0.1950 | 0.2010 | 0.192759 |
| 0.6 | 0.4530 | 0.4711 | 0.471513 |
| 0.7 | 0.7920 | 0.7834 | 0.774040 |
| 0.8 | 0.9590 | 0.9493 | 0.948097 |
| 0.9 | 0.9940 | 0.9965 | 0.995711 |
| 0.95 | 1.0000 | 0.9997 | 0.999558 |
| 0.99 | 1.0000 | 1.0000 | 0.999997 |

סעיף E (p=0.9, M=1000, חמש הרצות): ממוצע 0.9962, טעות יחסית 0.13% (תיאורטית 0.20% לפי האומדן).

עמודת "R מדויק" מגיעה מקובץ האימות ואינה מופיעה בדוח המוגש. כל אומדני ה-CMC עם M=10000 נמצאים במרחק של פחות מ-0.01 ממנה.

## חלק ב' – Destruction Spectrum מול Crude Monte Carlo

אותה רשת ואותם טרמינלים (**T1=4, T2=10, T3=17**). התוכנה מעריכה את ה-D-spectrum בעזרת תמורות אקראיות
(העוגן נמצא בעזרת DSS, בבנייה של הרשת מהסוף להתחלה), מחשבת ממנו את האמינות לכל p, ומשווה ל-CMC.

```
part-b-dspectrum/
├── assignment/project_part_b.pdf     # דף העבודה של חלק ב'
├── dspectrum_reliability.py          # הפתרון: D-spectrum + CMC, סעיפים (1)-(6)
├── validation/verify_part_b.py       # בדיקות אימות (לא חלק מההגשה)
├── results/
│   ├── run_output.txt                # טבלאות 1-4 (seed=2026)
│   └── verify_output.txt             # פלט בדיקות האימות
└── report/
    ├── אמינות_רשתות_פרויקט_חלק_ב.pdf  # הדוח להגשה
    ├── אמינות_רשתות_פרויקט_חלק_ב.docx
    ├── report_template.md, make_report.py  # הדוח נבנה מהפלט, בלי העתקה ידנית של מספרים
    └── build.py
```

הרצה (כ-8 שניות):

```bash
cd part-b-dspectrum
python dspectrum_reliability.py
python validation/verify_part_b.py      # אופציונלי
```

תוצאה מרכזית (טבלה 4, p=0.95, M=1000, עשר חזרות): טעות יחסית 1.74e-04 ב-D-spectrum מול 8.50e-04 ב-CMC.
