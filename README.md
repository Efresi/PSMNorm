
# OMOP_SEPSIS

Questo repository contiene il codice per condurre le analisi per lo studio "Predittori per la sepsi ospedaliera. Uno studio di coorte longitudinale del nodo OHDSI italiano.". 

Si prega di seguire le istruzioni sotto riportate per eseguire correttamente le analisi necessarie.

## Istruzioni

### Configurazione

Per eseguire le analisi è necessario avere installato R, RTools (Windows), RStudio e Java. Si consiglia di installare R 4.2.3 per una massima compatibilità.

Le istruzioni per poter configurare correttamente l'ambiente R si trovano al seguente link: https://ohdsi.github.io/Hades/rSetup.html .

### Download del materiale per eseguire le analisi

Dopo aver installato e verificato la corretta installazione dei tool richiesti, scaricare il contenuto di questo repository sulla macchina su cui saranno eseguite le analisi. Per farlo, cliccare sul pulsante verde "<> Code" che si trova nella parte superiore della pagina, da cui sarà possibile:
* scaricare un file zip contentente il materiale incluso nel repository, oppure
* clonare il repository utilizzando git.

Il repository contiene due cartelle principali:
* [OMOP_SEPSIS]-CohortDiagnostics
* [OMOP_SEPSIS]-PLPAnalysis

Ogni cartella contiene il codice per eseguire le analisi della corrispondente fase dello studio. All'interno di ogni cartella è presente un file README.md, che contiene le istruzioni dettagliate per eseguire correttamente il codice.

### Flusso di lavoro

Lo studio è diviso in due fasi principali:
1. esecuzione del tool CohortDiagnostics, per valutare le coorti incluse nello studio e la disponibilità dei dati nei centri partecipanti;
2. validazione del modello predittivo, sviluppato in precedenza dal centro coordinatore Fondazione IRCCS Policlinico San Matteo di Pavia, in collaborazione con l’Università di Pavia (UNIPV).

L'inizio di ognuna delle fasi è comunicato tramite e-mail dal centro coordinatore, in collaborazione con UNIPV.

Dopo la ricezione della prima e-mail, avrà inizio la fase 1 dello studio. Il codice che deve essere eseguito durante questa fase è quello che si trova all'interno della cartella [OMOP_SEPSIS]-CohortDiagnostics. 
Una volta che il codice è stato eseguito e i risultati condivisi, secondo le modalità comunicate dal coordinatore, si dovrà attendere la ricezione della seconda e-mail prima di procedere con la fase 2. 

Dopo aver ricevuto la seconda e-mail da parte del centro coordinatore, si potrà procedere con l'esecuzione del codice contenuto nella cartella [OMOP_SEPSIS]-PLPAnalysis.

## Problemi

Se riscontri problemi o hai domande, per favore apri una segnalazione qui su GitHub sfruttando la funzionalità **Issues** dedicata. Altri centri potrebbero avere lo stesso problema e in questo modo è possibile facilitarne la risoluzione per tutti.

In alternativa, invia una mail a e.fresi@smatteo.pv.it e/o a vittoria.ramella@biomeris.it.
