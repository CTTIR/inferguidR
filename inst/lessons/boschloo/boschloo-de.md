# Zwei unabhängige Binomialgruppen: Boschloo-Test

Diese eigenständige Lektion verwendet synthetische Anzahlen und das zugehörige Skript `boschloo.R`. Sie behandelt ausdrücklich die **zweiseitige Minlike-Variante**. Sie bestätigt keine vollständige numerische Freigabezertifizierung des Pakets.

## Fragestellung und Design

Zwei unabhängige Gruppen haben feste Gruppengrößen. Jede Person trägt genau eine binäre Beobachtung bei: Ereignis oder kein Ereignis. Unterscheiden sich die Ereigniswahrscheinlichkeiten? Ereignisdefinition, Gruppenreihenfolge, zweiseitige Alternative und Signifikanzniveau werden vor Sichtung der Anzahlen festgelegt. Gepaarte Beobachtungen, wiederholte Ereignisse derselben Person und komplexe Stichprobengewichte erfordern eine andere Analyse.

Gruppe 1 hat hier 7 Ereignisse bei 8 Personen, Gruppe 2 eines bei 8 Personen. Die beobachtete Risikodifferenz, Gruppe 1 minus Gruppe 2, beträgt 0,75. Das ist ein Schätzwert und keine Annahme für eine künftige Studienplanung.

## Test und Variante

Boschloo ordnet mögliche Tabellen mithilfe eines bedingten Fisher-p-Wertes. Anschließend wird die Nullwahrscheinlichkeit gleich oder stärker extremer Tabellen über die gemeinsame Ereigniswahrscheinlichkeit maximiert. Das entspricht dem Design zweier unabhängiger Binomialstichproben, ohne die Analyse auf die beobachtete Gesamtzahl der Ereignisse zu bedingen. Siehe die [Originalarbeit](https://doi.org/10.1111/j.1467-9574.1970.tb00104.x).

Die zweiseitigen Varianten **minlike** und **central** unterscheiden sich. Minlike verwendet die zweiseitige Fisher-Ordnung; central kombiniert die getrennt maximierten einseitigen Tests. Analyse und Planung dieser Lektion verwenden durchgehend minlike. Bei einer tatsächlich vorab festgelegten gerichteten Frage bedeutet `greater` p1 > p2 und `less` p1 < p2. Die Richtung darf nicht nach Sichtung des Ergebnisses geändert werden. Auch die [Dokumentation von exact2x2](https://search.r-project.org/CRAN/refmans/exact2x2/html/boschloo.html) unterscheidet die zweiseitigen Konventionen.

## Erst vorhersagen, dann ausführen

Müssen der unbedingte p-Wert und der bedingte Fisher-p-Wert gleich sein? Wie verändert sich die zukünftige Power, wenn p1 von 0,8 auf 0,7 sinkt und p2 bei 0,2 bleibt? Notieren Sie Ihre Antworten, bevor Sie die Ausgabe ansehen.

In R mit installiertem inferguidR:

```r
lesson <- new.env()
sys.source(inferguidR::learning_lesson("boschloo", "de")$script, envir = lesson)
answer <- lesson$run_boschloo_lesson(tempfile("boschloo-example-"))
answer
```

Das Skript protokolliert Paket-/R-Versionen und Sitzungsinformationen. Es schreibt ausschließlich in das angegebene Ausgabeverzeichnis. Englische/deutsche JSON- und TXT-Berichte werden erzeugt; PDF-Berichte entstehen bei verfügbarem Cairo, ohne automatische Installation eines Renderers. Hinzu kommen Planungs- und Sensitivitätstabellen sowie eine separate Planungsgrafik.

## Ergebnisse richtig lesen

Das ausgeführte Beispiel ergibt ungefähr Boschloo p=0,004180908 und Fisher p=0,01010101. Beide Tests verwerfen hier bei 0,05. Daraus folgt keine Austauschbarkeit zweiseitiger Konventionen. Vergleichen Sie Stichprobenannahmen und festgelegte Variante, nicht bloß den kleineren Zahlenwert.

`p_upper` begrenzt den Fehler der numerischen Optimierung. Es ist **keine** Stichproben-Konfidenzgrenze für p-Wert, Risikodifferenz oder Ereigniswahrscheinlichkeit. Das Skript berechnet getrennt zwei nicht gemeinsam adjustierte marginale 95%-Clopper–Pearson-Intervalle: ungefähr [0,47349; 0,99684] und [0,00316; 0,52651]. Diese bilden weder ein simultanes Intervallpaar noch ein Konfidenzintervall der Differenz. Ihre Überlappung testet nicht die Gleichheit der Wahrscheinlichkeiten. Ein Stichproben-Konfidenzintervall der Risikodifferenz wird in diesem Beispiel nicht berechnet.

## Künftige Studie mit derselben Variante planen

Die zukünftigen Annahmen sind p1=0,8, p2=0,2, alpha=0,05 und Power=0,8 bei gleichen Gruppengrößen. Diese hypothetischen Designwerte unterscheiden sich bewusst von den beobachteten Anteilen. Die begrenzte Suche reicht bis Gesamt-N=40 und prüft drei nachfolgende Designs auf lokale Stabilität.

Das erste passende Design hat n1=n2=10, insgesamt N=20, und ungefähr Power=0,8053902. Auch das festgelegte Stabilitätsfenster wird erfüllt. Das Minimum gilt entlang der angegebenen Allokationsfolge und im Suchbereich. Das lokale Stabilitätsfenster beweist nicht, dass jedes größere Design das Ziel erreicht. Die Power exakter Tests kann diskrete Sprünge haben. Deshalb zeigt die Grafik berechnete Designs und keine interpolierte Stichprobengrößengarantie.

Beim selben 10+10-Design beträgt die Power ungefähr 0,621632, 0,805390 und 0,936660 für p1=0,7, 0,8 und 0,9 bei festem p2=0,2. Das sind Planungsszenarien, keine Konfidenzintervalle der beobachteten Daten und kein rückblickender Beleg für die Hypothese.

## Verwandte Lehrmaterialien und ausdrücklich offene Punkte

Das [CTTIR-Fisher-Tutorial](https://cttir.github.io/tutorials/tutorials/inference/fisher-exact-test.html) dient **nur als verwandter Hintergrund**, nicht als Boschloo-Lektion. Seine angezeigten Zahlen passen nicht zur dargestellten Tabelle (bekanntes Problem SP10). Führen Sie den Code erneut aus, statt die Zahlen zu übernehmen. [Festgehaltene Quellversion als Ersatzlink](https://github.com/CTTIR/tutorials/blob/7d8d01c4edabd294f88568c4e370ab83dc09d313/tutorials/inference/fisher-exact-test.qmd). Im geprüften Quellbestand wurde keine CTTIR-Boschloo-Lektion gefunden. Diese eigenständige lokale Lektion schließt genau diese Lücke; weitere Lehrthemen bleiben gesonderte Aufgaben.

Ein Bericht soll Stichprobendesign, Ereignisrichtung, Variante, beobachteten Effekt, p-Wert, numerischen Status und fehlendes Intervall ausdrücklich nennen. Ein erfolgreicher Lehrdurchlauf ist keine vollständige wissenschaftliche Zertifizierung.
