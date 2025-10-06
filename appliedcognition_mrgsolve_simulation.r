# Copyright (c) 2025 Applied Cognition Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# mrgSolve glymphatic flow compartment simulation

rm(list = ls())
library(mrgsolve)

file_dir = "./output"
dir.create(file_dir, showWarnings = F)

# ---- Model ----

# Model notes
# - PRODUCTION is implemented via EV infusion (amt/rate) into cmt 1 (CELL).
# - CLEARANCE scales K12 and K21 (glymphatic-like effect).
# - Ab42/pTau enable oligomeric dynamics; Ab40/nTau do not.

model_code <- '
$PARAM
CLEARANCE     = 1
PRODUCTION    = 1
OLIGO_ON_RATE = 0
OLIGO_OFF_RATE= 0
Ab40          = 0
nTau          = 0
pTau          = 0

$INIT
CELL        = 0
ISF         = 0
CSF         = 0
PLASMA      = 0
NEUROPIL    = 0
OLIGOMER    = 0

$MAIN
// Base rate constants
double K120 = 0.151;
double K210 = 0.151;
double K13  = 0.1*0.1;
double K14  = 0.1*0.9;
double K23  = 0.151;
double K3e  = 1.0;
double K4e  = 1.0;

// Sleep/wake Production dynamics
double K01  = PRODUCTION;

// Oligomer dynamics
double K15  = OLIGO_ON_RATE;
double K51  = OLIGO_OFF_RATE;

// Apply CLEARANCE scaling
double K12 = K120 * CLEARANCE;
double K21 = K210 * CLEARANCE;

// Tau shuts off BBB transport
if ((nTau == 1) || (pTau == 1)) {
  K13 = 0;
}

$ODE
dxdt_CELL        = -K01 * (CELL-100);
dxdt_ISF         = -(K15 + K12 + K13 + K14) * ISF + K51 * OLIGOMER + K21 * CSF + K01 * CELL;
dxdt_CSF         = -(K21 + K23) * CSF + K12 * ISF;
dxdt_PLASMA      = -K3e * PLASMA + K13 * ISF + K23 * CSF;
dxdt_NEUROPIL    = -K4e * NEUROPIL + K14 * ISF;
dxdt_OLIGOMER    = -K51 * OLIGOMER + K15 * ISF;

$TABLE
capture K12_t  = K12;
capture K21_t  = K21;
capture CLEAR_ = CLEARANCE;
capture PROD_  = PRODUCTION;
'

# Save and compile
dir.create(file.path(file_dir, "model"), showWarnings = FALSE, recursive = TRUE)
writeLines(model_code, file.path(file_dir, "model", "biomarker_comp_model.cpp"))
mod <- mread("biomarker_comp_model", file.path(file_dir, "model"))

# ---- Setup grid/time ----
time <- seq(0, 2400, by = 4)
Tlen <- length(time)

base <- data.frame(
  ID   = 1L,
  time = time,
  cmt  = 1L,
  ss   = 0L,
  evid = 1L
)

# ---- Condition coding (use simple named lists, not hash) ----
vc <- list(
  "Physio+" = list(clears = c(0.6, 0.8, 1, 1.2, 1.4), prods = c(1, 1, 1, 1, 1)),
  "Physio-" = list(clears = c(0.9, 0.95, 1, 1.05, 1.1), prods = c(1, 1, 1, 1, 1)),
  "Hypno+"  = list(clears = c(0.6, 0.8, 1, 1.2, 1.4), prods = c(1, 1, 1, 1, 1)),
  "Hypno-"  = list(clears = c(0.95, 0.975, 1, 1.025, 1.05), prods = c(1, 1, 1, 1, 1)),
  "Wake+"   = list(clears = rep(1, 5), prods = c(0.8, 0.9, 1, 1.1, 1.2)),
  "Wake-"   = list(clears = rep(1, 5), prods = c(0.8, 0.9, 1, 1.1, 1.2)),
  "Increase"= list(clears = c(0.67, 0.67, 1.00), prods = c(1, 1.25, 1)),
  "Decrease"= list(clears = c(0.67, 0.67, 0.33), prods = c(1, 0.5, 1)),
  "Extreme" = list(clears = c(1, 1, 1, 0.5, 2), prods = c(1.5, 0.5, 1, 1, 1))
)

biomks        <- c("Ab40","Ab42","nTau","pTau")
biomks_label  <- c('CM Ab40', 'CM Ab42','CM nTau', 'CM pTau')
conds         <- names(vc)
Models        <- c("null_model", "glymphatic")

# Labels are only used for the 5-step conditions in your overnight-change tables
change_labels5 <- c('-1','-0.5','0','0.5','1')

bio_setup <- function(bio) {
  #' Sets bio-specific flags & on/off rates

  out <- data.frame(
    Ab40 = as.integer(bio == "Ab40"),
    nTau = as.integer(bio == "nTau"),
    pTau = as.integer(bio == "pTau")
  )
  if (bio %in% c("Ab42", "pTau")) {
    out$OLIGO_ON_RATE  <- 7.2
    out$OLIGO_OFF_RATE <- 0.32
  } else {
    out$OLIGO_ON_RATE  <- 0
    out$OLIGO_OFF_RATE <- 0
  }
  out[rep(1, Tlen), , drop = FALSE]
}

make_series <- function(hypothesis, cond, clear, prod) {
  #' Creates CLEAR/PROD time series by hypothesis/condition/values

  # 6 equal blocks across the day; sleep = last two blocks
  rep6 <- function(v) rep(v, (Tlen - 1L) / 6L)
  if (hypothesis == "glymphatic") {
    if (cond %in% c("Wake+","Wake-")) {
      CLEAR <- 0.25 * c(1, rep6(c(1,1,1,1,1.1*clear,1.1*clear)))
    } else {
      CLEAR <- 0.25 * c(1, rep6(c(1,1,1,1,1.5*clear,1.5*clear)))
    }
  } else {
    CLEAR <- 0.25 * c(1, rep6(c(1,1,1,1,1,1)))
  }

  # Baseline ~30% drop during sleep; Wake+/- uses 0.9*prod else 0.7*prod
  if (cond %in% c("Wake+","Wake-")) {
    PROD <- 0.1 * c(1, rep6(c(1,1,1,1,0.9*prod,0.9*prod)))
  } else {
    PROD <- 0.1 * c(1, rep6(c(1,1,1,1,0.7*prod,0.7*prod)))
  }
  list(CLEAR = CLEAR, PROD = PROD)
}

# ---- Create Containers ----
plot_models <- list()

Hypothesis <- Condition <- Biomarker <- Biomarker_label <- character()
Clear <- Prod <- Change_label <- character()
Plasma_Post <- CSF_Post <- numeric()
Plasma_Overnight_Change <- CSF_Overnight_Change <- numeric()

# ---- Main loops (kept, but much leaner inside) ----
for (hypothesis in Models) {
  for (cond in conds) {
    clears <- vc[[cond]]$clears
    prods  <- vc[[cond]]$prods

    for (bio in biomks) {
      bio_block <- bio_setup(bio)

      for (k in seq_along(clears)) {
        clear <- clears[k];  prod <- prods[k]
        ser   <- make_series(hypothesis, cond, clear, prod)

        physio <- cbind(
          base, bio_block,
          amt        = 10 * c(rep(rep(4, 6), (Tlen - 1L) / 6L), 4),
          rate       = 10 * c(rep(rep(1, 6), (Tlen - 1L) / 6L), 1),
          CLEARANCE  = ser$CLEAR,
          PRODUCTION = ser$PROD
        )

        out <- mrgsim(mod, data = physio, end = 2400, delta = 0.1, maxsteps = 100000)

        if (cond %in% c("Increase","Decrease")) {
          if (hypothesis == "glymphatic") {
            akey <- sprintf("%s %s clear:%s prod:%s", cond, bio, clear, prod)
            plot_models[[akey]] <- out
          }
        } else {
          # Overnight change (use ss = 97 )
          ss <- 97
          t_pre  <- which(out$time == (24*ss + 16))[2]
          t_post <- which(out$time == (24*ss + 24))[2]

          Plasma_pre  <- out$PLASMA[t_pre]
          Plasma_post <- out$PLASMA[t_post]
          CSF_pre     <- out$CSF[t_pre]
          CSF_post    <- out$CSF[t_post]

          Hypothesis               <- c(Hypothesis, hypothesis)
          Condition                <- c(Condition,  cond)
          Biomarker                <- c(Biomarker,  bio)
          Biomarker_label          <- c(Biomarker_label, biomks_label[match(bio, biomks)])
          Clear                    <- c(Clear,  clear)
          Prod                     <- c(Prod,   prod)
          Change_label             <- c(Change_label, change_labels5[pmin(k, length(change_labels5))])
          Plasma_Post              <- c(Plasma_Post,  Plasma_post)
          CSF_Post                 <- c(CSF_Post,     CSF_post)
          Plasma_Overnight_Change  <- c(Plasma_Overnight_Change, (Plasma_post - Plasma_pre) / Plasma_pre)
          CSF_Overnight_Change     <- c(CSF_Overnight_Change,    (CSF_post    - CSF_pre)    / CSF_pre)
        }
      }
    }
  }
}

df_compart <- data.frame(
  Hypothesis, Condition, Biomarker, Biomarker_label,
  Clear = as.numeric(Clear), Prod = as.numeric(Prod),
  Change_label,
  Plasma_Post, CSF_Post,
  Plasma_OV_Change = Plasma_Overnight_Change,
  CSF_OV_Change    = CSF_Overnight_Change,
  row.names = NULL, stringsAsFactors = FALSE
)

# ---- Percent change vs label==0 within each (Hypothesis, Condition, Biomarker) ----
PrctChg_Plasma_Post <- numeric(0)
PrctChg_CSF_Post    <- numeric(0)

for (hypo in unique(df_compart$Hypothesis)) {
  for (cond in setdiff(conds, c("Increase","Decrease"))) {
    idx1 <- df_compart$Hypothesis == hypo & df_compart$Condition == cond
    df_  <- df_compart[idx1, ]
    for (b in unique(df_$Biomarker)) {
      idx2 <- df_$Biomarker == b
      base_row <- which(df_$Change_label[idx2] == '0')
      # guard: if missing baseline, skip
      if (length(base_row) == 0L) {
        PrctChg_Plasma_Post <- c(PrctChg_Plasma_Post, rep(NA_real_, sum(idx2)))
        PrctChg_CSF_Post    <- c(PrctChg_CSF_Post,    rep(NA_real_, sum(idx2)))
      } else {
        normP <- df_$Plasma_Post[idx2][base_row[1]]
        normC <- abs(df_$CSF_Post[idx2][base_row[1]])
        valsP <- df_$Plasma_Post[idx2]
        valsC <- df_$CSF_Post[idx2]
        PrctChg_Plasma_Post <- c(PrctChg_Plasma_Post, 100 * (valsP - normP) / abs(normP))
        PrctChg_CSF_Post    <- c(PrctChg_CSF_Post,    100 * (valsC - normC) / abs(normC))
      }
    }
  }
}

df_compart_ <- cbind(df_compart, PrctChg_Plasma_Post, PrctChg_CSF_Post)

# ---- Ratio changes for Extreme condition (kept close to original but simplified) ----
mk_ratio_block <- function(df_slice) {
  get_ref <- function(mk) with(df_slice, Clear == 1 & Prod == 1 & Biomarker == mk)
  get_mk  <- function(mk) df_slice$Biomarker == mk

  ref_Ab42P <- df_slice$Plasma_Post[get_ref("Ab42")]
  ref_Ab40P <- df_slice$Plasma_Post[get_ref("Ab40")]
  ref_Ab42C <- df_slice$CSF_Post[get_ref("Ab42")]
  ref_Ab40C <- df_slice$CSF_Post[get_ref("Ab40")]

  ref_pTauP <- df_slice$Plasma_Post[get_ref("pTau")]
  ref_nTauP <- df_slice$Plasma_Post[get_ref("nTau")]
  ref_pTauC <- df_slice$CSF_Post[get_ref("pTau")]
  ref_nTauC <- df_slice$CSF_Post[get_ref("nTau")]

  norm_Ab_P <- ref_Ab42P / ref_Ab40P
  norm_Ab_C <- ref_Ab42C / ref_Ab40C
  norm_T_P  <- ref_pTauP / ref_nTauP
  norm_T_C  <- ref_pTauC / ref_nTauC

  i40 <- get_mk("Ab40")  # use Ab40 to carry Clear/Prod lengths
  data.frame(
    Hypothesis = df_slice$Hypothesis[i40],
    Condition  = "Extreme",
    Clear      = df_slice$Clear[i40],
    Prod       = df_slice$Prod[i40],
    Plasma_Ratio.Ab42.Ab40 = (df_slice$Plasma_Post[get_mk("Ab42")] / df_slice$Plasma_Post[get_mk("Ab40")]) / as.numeric(norm_Ab_P),
    CSF_Ratio.Ab42.Ab40    = (df_slice$CSF_Post[get_mk("Ab42")]    / df_slice$CSF_Post[get_mk("Ab40")])    / as.numeric(norm_Ab_C),
    Plasma_Ratio.pTau.nTau = (df_slice$Plasma_Post[get_mk("pTau")] / df_slice$Plasma_Post[get_mk("nTau")]) / as.numeric(norm_T_P),
    CSF_Ratio.pTau.nTau    = (df_slice$CSF_Post[get_mk("pTau")]    / df_slice$CSF_Post[get_mk("nTau")])    / as.numeric(norm_T_C)
  )
}

df_ratios <- do.call(rbind, lapply(Models, function(hypo) {
  mk_ratio_block(subset(df_compart, Condition == "Extreme" & Hypothesis == hypo))
}))

write.csv( df_compart_, file = file.path(file_dir, "df_compart_BIO_P2_SIMULATION.csv"), row.names = TRUE)
write.csv( df_ratios, file = file.path(file_dir, "df_compart_BIO_P2_SIMULATION_RATIOS.csv"), row.names = TRUE)



# ---- Plots ----

plot_panels <- function(bio, increase = TRUE, output_format = 'pdf') {
  #' Generates the three-panel plot: Plasma (top), CSF (middle), ISF & Cell (bottom)

  if (increase) {
    fig          <- c(paste("Increase", bio,"clear:0.67 prod:1"),
                      paste("Increase", bio, "clear:1 prod:1"),
                      paste("Increase", bio, "clear:0.67 prod:1.25"))
    dotted_line  <- "Increased Production"
    dashed_line  <- "Increased Clearance"
    fig_title = paste0("INCREASE_", bio)
  } else {
    fig          <- c(paste("Decrease", bio, "clear:0.67 prod:1"),
                      paste("Decrease", bio, "clear:0.33 prod:1"),
                      paste("Decrease", bio, "clear:0.67 prod:0.5"))
    dotted_line  <- "Decreased Production"
    dashed_line  <- "Decreased Clearance"
    fig_title = paste0("DECREASE_", bio)
  }

  # Grab runs
  m_null <- plot_models[[fig[1]]]
  m_clr  <- plot_models[[fig[2]]]
  m_prd  <- plot_models[[fig[3]]]

  # Window: last 48 h (two 8 h "sleep" windows as in your code)
  t0    <- 2400 - 72
  tend  <- t0 + 48
  idx   <- function(m) (m$time >= t0) & (m$time <= tend)
  ix0   <- idx(m_null)
  ixc   <- idx(m_clr)
  ixp   <- idx(m_prd)

  # Consistent x ticks labeled as hours since t0
  xticks   <- seq(t0, tend, by = 8)
  xlabels  <- as.character(seq(0, 48, by = 8))

  # Sleep rectangles (two 8-hour windows)
  sleep_windows <- list(
    c(t0 + 16, t0 + 24),
    c(t0 + 40, t0 + 48)
  )

  shade_sleep <- function() {
    usr <- par("usr")
    for (sw in sleep_windows) {
      rect(xleft = sw[1], xright = sw[2], ybottom = usr[3], ytop = usr[4],
           col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
    }
  }

  add_grid <- function() {
    op <- par(xpd = FALSE)
    on.exit(par(op), add = TRUE)
    abline(v = xticks, col = "darkgray", lty = "dotted")
    abline(h = axTicks(2), col = "darkgray", lty = "dotted")
  }

  # Line types per scenario
  lty_null <- 1
  lty_clr  <- 2
  lty_prd  <- 3

  # Colors per compartment
  col_plasma <- "red"
  col_csf    <- "blue"
  col_isf    <- "darkgreen"
  col_cell   <- "magenta"

  # Compute y-lims per panel across the three scenarios within the window
  yrng_panel <- function(values_list) {
    rng <- range(unlist(values_list), finite = TRUE)
    if (!all(is.finite(rng))) rng <- c(0, 1)
    if (diff(rng) == 0) rng <- rng + c(-0.5, 0.5)
    rng
  }

  yl_plasma <- yrng_panel(list(m_null$PLASMA[ix0], m_clr$PLASMA[ixc], m_prd$PLASMA[ixp]))
  yl_csf    <- yrng_panel(list(m_null$CSF[ix0],    m_clr$CSF[ixc],    m_prd$CSF[ixp]))
  yl_rest   <- yrng_panel(list(m_null$ISF[ix0],    m_clr$ISF[ixc],    m_prd$ISF[ixp],
                               m_null$CELL[ix0],   m_clr$CELL[ixc],   m_prd$CELL[ixp]))

  # --- sizing controls ---
  AXIS_LAB_CEX <- 1.5*0.8   # y/x axis labels
  TICK_CEX     <- 1.4 *0.8 # tick labels
  LEGEND_CEX   <- 1.4*0.8  # legend text
  MAIN_CEX     <- 1.25*0.8  # main/outer titles (if any)

  # set global sizes
  par(cex.lab = AXIS_LAB_CEX,
      cex.axis = TICK_CEX,
      cex.main = MAIN_CEX)
  
  if (output_format == 'pdf') {
    # Save to PDF
    pdf(file = file.path(file_dir, paste0(fig_title, '.pdf')),
        width = 6, height = 6, pointsize = 10,
        onefile = FALSE,        # multiple pages if you make multiple plots
        useDingbats = FALSE)
  } else if (output_format == 'tiff') {
    # Save to TIFF
    tiff(file = file.path(file_dir, paste0(fig_title, '.tiff')),
         width = 6, height = 6, pointsize = 10, units='in', res=300)
  } else {
    print('Unknown plot output format')
    quit(status=1)
  }

  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)

  # lay out 3 horizontal panels on ONE page
  par(mfrow = c(3, 1), mar = c(2.5, 5, 2.5, 1.2), oma = c(2, 0, 1.8, 0))

  # ---- Panel 1: PLASMA ----
  plot(m_null$time[ix0], m_null$PLASMA[ix0], type = "n",
       xlab = "", xaxt = "n",
       ylab = "Plasma [a.u.]", ylim = yl_plasma)

  shade_sleep()
  add_grid()

  lines(m_null$time[ix0], m_null$PLASMA[ix0], col = col_plasma, lwd = 2, lty = lty_null)
  lines(m_clr$time[ixc],  m_clr$PLASMA[ixc],  col = col_plasma, lwd = 2, lty = lty_clr)
  lines(m_prd$time[ixp],  m_prd$PLASMA[ixp],  col = col_plasma, lwd = 2, lty = lty_prd)

  axis(1, at = xticks, labels = xlabels, cex.axis = TICK_CEX)
  #mtext("Time [hr]", side = 1, line = 2, cex = AXIS_LAB_CEX)
  legend("topleft", inset = 0.01,
         legend = c("Plasma"),
         col    = c(col_plasma),
         lty    = c(1), lwd = 2, bty = "n", cex = LEGEND_CEX)


  # Global legend for scenario line types (null / clearance / production)
  par(xpd = NA)
  legend("top", inset = -0.3, horiz = TRUE,
         legend = c("Null (solid)", dashed_line, dotted_line),
         lty    = c(lty_null, lty_clr, lty_prd),
         col    = "darkgray", lwd = 2, bty = "n", cex = 1.2*LEGEND_CEX)
  par(xpd = FALSE)  # restore clipping

  # ---- Panel 2: CSF ----
  plot(m_null$time[ix0], m_null$CSF[ix0], type = "n",
       xlab = "", xaxt = "n",
       ylab = "CSF [a.u.]", ylim = yl_csf)

  shade_sleep()
  add_grid()

  lines(m_null$time[ix0], m_null$CSF[ix0], col = col_csf, lwd = 2, lty = lty_null)
  lines(m_clr$time[ixc],  m_clr$CSF[ixc],  col = col_csf, lwd = 2, lty = lty_clr)
  lines(m_prd$time[ixp],  m_prd$CSF[ixp],  col = col_csf, lwd = 2, lty = lty_prd)

  axis(1, at = xticks, labels = xlabels, cex.axis = TICK_CEX)
  #mtext("Time [hr]", side = 1, line = 2, cex = AXIS_LAB_CEX)
  legend("topleft", inset = 0.01,
         legend = c("CSF"),
         col    = c(col_csf),
         lty    = c(1), lwd = 2, bty = "n", cex = LEGEND_CEX)

  # ---- Panel 3: ISF & CELL ----
  plot(m_null$time[ix0], m_null$ISF[ix0], type = "n",
       xlab = "", xaxt = "n",
       ylab = "Other [a.u.]", ylim = yl_rest)

  shade_sleep()
  add_grid()

  # ISF
  lines(m_null$time[ix0], m_null$ISF[ix0], col = col_isf,  lwd = 2, lty = lty_null)
  lines(m_clr$time[ixc],  m_clr$ISF[ixc],  col = col_isf,  lwd = 2, lty = lty_clr)
  lines(m_prd$time[ixp],  m_prd$ISF[ixp],  col = col_isf,  lwd = 2, lty = lty_prd)
  # CELL
  lines(m_null$time[ix0], m_null$CELL[ix0], col = col_cell, lwd = 2, lty = lty_null)
  lines(m_clr$time[ixc],  m_clr$CELL[ixc],  col = col_cell, lwd = 2, lty = lty_clr)
  lines(m_prd$time[ixp],  m_prd$CELL[ixp],  col = col_cell, lwd = 2, lty = lty_prd)
  axis(1, at = xticks, labels = xlabels, cex.axis = TICK_CEX)
  mtext("Time [hr]", side = 1, line = 3, cex = AXIS_LAB_CEX)

  # Panel legend: which compartments are shown here
  legend("topleft", inset = 0.01,
         legend = c("ISF", "Cell"),
         col    = c(col_isf, col_cell),
         lty    = c(1, 1), lwd = 2, bty = "n", cex = LEGEND_CEX)

  dev.off()

  return(0)
}


# ---- Create plots  ----

bio = "Ab40"
plot_panels(bio, increase = TRUE, output_format = 'pdf')
plot_panels(bio, increase = FALSE, output_format = 'pdf')

plot_panels(bio, increase = TRUE, output_format = 'tiff')
plot_panels(bio, increase = FALSE, output_format = 'tiff')
